%function [res_gen, res_pol] = resultPS(mpc, esc, idx_gen, gen_wgt, cur_info, pol_val, pol_name, res, res_gen, res_pol)

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%%
idx_stor = strcmp(mpc.genfuel, 'storage');
% if any indices in idx_gen refers to demand side
if any(idx_gen & strcmp(mpc.gentype, 'dl'))
    % Remove things purely on demand side from idx_gen
    idx_gen = idx_gen & ~(ismember(mpc.genfuel, {'dl', 'storage'}));
    idx_dl = gen_wgt .* ismember(mpc.genfuel, {'dl', 'dac', 'storage'});
    pol_val = sum(gen_wgt(idx_dl > 0).*mpc.gen(idx_dl > 0, 2)) / sum(mpc.gen(idx_dl > 0, 2));
else
    idx_area_gen = getInfoIdx(mpc, esc, rmfield(cur_info, {'genfuel', 'gentype'}), 'gen', 1);
    idx_dl = pol_val * ones(mpc.ng, 1) .* idx_area_gen .* ismember(mpc.genfuel, {'dl', 'dac', 'storage'});
    idx_dl(idx_stor) = idx_dl(idx_stor) * -0.15/0.85;
    idx_gen = idx_gen | idx_area_gen .* strcmp(mpc.genfuel, 'dac');
    %pol_val = 1;
end
if ~(pol_val >= 0 && pol_val <= 1) % target value, convert to fraction
    idx_area_bus = getInfoIdx(mpc, esc, rmfield(cur_info, {'genfuel', 'gentype'}), 'bus');
    pol_val = pol_val / sum(mpc.ann_load(idx_area_bus));
end

idx_mu = strcmp(mpc.total_output.name, char(pol_name));
% If the total output constraint doesn't get made for any reason (all zeros?), then the rest of the results will not work
if ~any(idx_mu)
    vfprintf(eopt.verbose, "Total output constraint for policy %s was not created. Skipping results for it \n", char(pol_name))
    return
end
col = mpc.total_output.type(idx_mu, :);

% in an rps, dac pays for their credits (in dl not gen)
% and in ces, dac should get paid for their credits (in gen not dl)
% but want to allow for dacs in ces not generating enough credits to cover
% its elec use, and has to pay for credits
if contains(pol_type, 'rps')
    % taking dac out of rps gens, for same reason as in applyPS.
    idx_gen = idx_gen & ~strcmp(mpc.genfuel, 'dac');
elseif contains(pol_type, 'ces')
    coeff = mpc.total_output.coeff(:, col);
    % remove dacs generating credits from dl
    idx_dl = idx_dl .* ~(strcmp(mpc.genfuel, 'dac') & coeff >= 0);
    % remove dacs consuming credits from gen
    idx_gen = idx_gen & ~(strcmp(mpc.genfuel, 'dac') & coeff <=0);
end

if ~(sum(idx_dl) == 0 || isnan(sum(pol_val)))

    map = abs(mpc.total_output.coeff(:, col)); % Coeff Sign is Flipped
    pol_prc = res.total_output.mu(idx_mu);
    % Right now, using generation_mwh for storage, which represents loss, so
    % revert the coefficients to be per-loss
    map(idx_stor) = map(idx_stor)/(0.15/0.85);
    
    res_gen{idx_dl ~= 0, [char(pol_name), '_req']} = idx_dl(idx_dl ~= 0);
    res_gen{map ~= 0, [char(pol_name), '_credits']} = abs(res_gen{map ~= 0, 'generation_mwh'}.*map(map ~= 0));
    % % For storage, generation_mwh is measure of actual loss, but TOC coefficient is in terms of storage's positive generation.
    %res_gen{map ~= 0 & idx_stor, [char(pol_name), '_credits']} = abs(res_gen{map ~= 0 & idx_stor, 'generation_pos_mwh'}.*map(map ~= 0 & idx_stor));
    res_gen{map ~= 0, [char(pol_name), '_mu']} = pol_prc;
    
    pol_qty = sum(res_gen{idx_gen, [char(pol_name), '_credits']});
    
    % When a ces/rps constraint isn't binding, dl and gens have different #
    % of credits, and policy isn't active, so we dont need to add it to 
    % welfare, etc. Usually an inactive constraint has mu < 0.0001 
    % i think we want no overlap between idx_gen and idx_dl > 0 
    %if abs(pol_prc) > 0.00001
    supply = sum(res_gen{idx_gen, [char(pol_name), '_credits']}); % = pol_qty
    demand = sum(res_gen{idx_dl ~= 0, [char(pol_name), '_credits']});
    surplus = supply - demand;
    
    if abs(surplus) < 25 || abs(pol_prc) > 0.00001
        res_gen{map ~= 0, [char(pol_name), '_cost']} = res_gen{map ~= 0, [char(pol_name), '_mu']} .* res_gen{map ~= 0, [char(pol_name), '_credits']};
        res_gen{idx_dl ~= 0, f_wel} = res_gen{idx_dl ~= 0, f_wel} - res_gen{idx_dl ~= 0, [char(pol_name), '_cost']};
        res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};
        
        res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_prc, (pol_qty * pol_prc)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];
        if abs(surplus) > 25
            vfprintf(1, "Warning: %s is binding, but %.3f credits purchased by units are unaccounted for. This is probably due to storage wasting.\n", pol_name{:}, surplus);
        end
    else
        res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_prc, 0], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];
    end
end
