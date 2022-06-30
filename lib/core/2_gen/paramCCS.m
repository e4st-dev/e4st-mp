function [mpc, offer] = paramCCS(mpc, esc, offer, eopt)
% paramCCS Set up CCS generators

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% NOTES
% Assumes 90% carbon capture for all retrofit units. Also assumes that the
% NOX emissions rates are cut in half when the plant is retrofit and SO2 
% emissions rates are zero. Both of these assumptions are made directly in the code.

%% Important parameters
define_constants;
dol_conv = 1 / 1.04107906245393; % 2016$ to 2013$
%capital recovery factor
if ~isfield(eopt, 'ccsretro_crf')
    eopt.ccsretro_crf = 0.0533;
end
% Storage options per plant
stor_steps = esc.stor_steps;
n_stor = eopt.n_stor;

%% Units to Retrofit
if strcmp(eopt.ccs_retro, 'T')
    idx_online = mpc.gen(:, GEN_STATUS) == 1;
    idx_retro = idx_online & ...
        (strcmp(mpc.gentype, 'coal') | strcmp(mpc.gentype, 'ngcc')) & ...
        strcmp(mpc.gen_map{:, 'nation'}, 'us');

    idx_retro = idx_retro & ~strcmp(mpc.gentype, 'ngcc'); % turn off ngcc retrofits
    
    %% Coal CCS Retrofits only if avgcap above 300 MW. 
    % This reduces the problem size and helps us avoid extrapolating the heat rate data
%     idx_retro = idx_retro & mpc.gen_desc{:, 'AVG_CAP'} >= 300;

    %% Initialize Unit Data
    [retrogen, retrogen.offer] = filter_e4st_gen(mpc, offer, idx_retro);
    
    
    %% Apply CCS Parameters
    
    retrogen.newgen(:, :) = 1;
    avgcap = retrogen.gen_desc{:, 'AVG_CAP'};
    hr = retrogen.gen_aux{:, 'HR'};
    
    % Coal (90% CCS):
    % Modeled after Petra Nova CCS
    % Steam for CCS system comes from auxillary NG unit: Small capacity
    % penalty, reduction in CCS capture; no heat rate penalty.
    idx_retro_coal = strcmp(retrogen.genfuel, 'coal') &  ~strcmp(retrogen.gentype, 'coalccs_new');
    retrogen.gentype(idx_retro_coal) = {'coalccs_ret'};  
    
    retrogen.gen_aux{idx_retro_coal, 'CAP_COST'} = retrogen.gen_aux{idx_retro_coal, 'CAP_COST'} + (ccs_retrofit_curves('CAP_COST', avgcap, hr) * dol_conv * eopt.ccsretro_crf / 8.76);
    retrogen.gen_aux{idx_retro_coal, 'FOM'} = retrogen.gen_aux{idx_retro_coal, 'FOM'} + (ccs_retrofit_curves('FOM', avgcap, hr) * dol_conv / 8.76);
    retrogen.gen_aux{idx_retro_coal, 'VOM'} = retrogen.gen_aux{idx_retro_coal, 'VOM'} + (ccs_retrofit_curves('VOM', avgcap, hr) * dol_conv);
    retrogen.gen_aux{idx_retro_coal, 'HR'} = retrogen.gen_aux{idx_retro_coal, 'HR'} .* (1 + (ccs_retrofit_curves('HRPen', avgcap, hr) / 100));
    retrogen.gen_aux{idx_retro_coal, 'EMIS_CO2'} = retrogen.gen_aux{idx_retro_coal, 'EMIS_CO2'} * (1 - .9) .* (1 + (ccs_retrofit_curves('HRPen', avgcap, hr) / 100)); % Apply HR penalty
    retrogen.gen_aux{idx_retro_coal, 'EMIS_NOX'} = retrogen.gen_aux{idx_retro_coal, 'EMIS_NOX'} / 2 .* (1 + (ccs_retrofit_curves('HRPen', avgcap, hr) / 100)); % Cut NOX in half Apply HR penalty;
    retrogen.gen_aux{idx_retro_coal, 'EMIS_SO2'} = 0;
    retrogen.gen_aux{idx_retro_coal, 'EMIS_PM25'} = retrogen.gen_aux{idx_retro_coal, 'EMIS_PM25'} * 0.65 .* (1 + (ccs_retrofit_curves('HRPen', avgcap, hr) / 100)); % Apply HR penalty;
    retrogen.gen(idx_retro_coal, PMAX) = retrogen.gen(idx_retro_coal, PMAX) .* (1 - (ccs_retrofit_curves('CapPen', avgcap, hr) / 100));
    retrogen.offer(idx_retro_coal, 2) = retrogen.offer(idx_retro_coal, 2) .* (1 - (ccs_retrofit_curves('CapPen', avgcap, hr) / 100));
    

    %% Add Retrofits to Main Data
    [mpc, offer] = append_e4st_gen(mpc, offer, retrogen);
    
    %% Update Offer and Gencost
    idx_retro = strcmp(mpc.gentype, 'coalccs_ret') | strcmp(mpc.gentype, 'ngccccs_ret');
    offer = updateOfferPrc(mpc, offer, idx_retro, eopt);
    mpc = updateGenCost(mpc, idx_retro, eopt);
    
    %% Calculate CO2 Storage Rate
    mpc.gen_aux{idx_retro, 'CO2_CAPTURE_PC'} = 0.9;
    mpc.gen_aux{idx_retro, 'STORAGE_CO2'} = mpc.gen_aux{idx_retro, 'EMIS_CO2'} ./ (1 - mpc.gen_aux{idx_retro, 'CO2_CAPTURE_PC'}) .* mpc.gen_aux{idx_retro, 'CO2_CAPTURE_PC'}; % 90% CCS, how much CO2 is stored
  
end


%% CCS Storage Options

idx_ccs = ismember(mpc.gentype, {'coalccs_ret', 'ngccccs_ret', 'ngccccs_new', 'coalccs_new'});
% We don't need to re-create storage options for units that already have
% them from previous iteration or anything else. Assuming default is STEP0
% Note that this doesn't allow CCS units to change their storage location
% between years, as existing ccs don't receive new options.
idx_ccs = idx_ccs & strcmp(mpc.gen_map{:, 'co2_storage_step'}, 'STEP0'); %only the unassigned CCS


gen_states = mpc.gen_map{:, 'state'};
prod_states = unique(gen_states(idx_ccs));
for prod_state = prod_states'
    idx_gen = idx_ccs & strcmp(gen_states, prod_state); %gens in the state
    gen_stor_steps = stor_steps(strcmp(stor_steps{:, 'ProdState'}, prod_state), :);
    gen_stor_steps = sortrows(gen_stor_steps, 'total_cost_2013'); %order from cheapest to most expensive
    gen_stor_steps = gen_stor_steps(1:n_stor, :); %only keep the first n
    if isfield(eopt, 'co2_stor_type')
        % If need to only have a certain storage type
        idx_keep = ismember(gen_stor_steps{:, 'stor_type'}, eopt.co2_stor_type);
        gen_stor_steps = gen_stor_steps(idx_keep, :);
    end
    [mpc_subset, mpc_subset.offer] = filter_e4st_gen(mpc, offer, idx_gen);
    stor_states = gen_stor_steps{:, 'StorState'}; % storage state options for the gens in prod_state
    fprintf('Creating CCS Storage Options by Production State: %s (%d generators x %d storage states = %d total)\n', ...
        char(prod_state), sum(idx_gen), length(stor_states), sum(idx_gen)*length(stor_states))
    for i_stor_state = 1:length(stor_states)
        stor_state = stor_states(i_stor_state);
        %fprintf('\tStorage State: %s\n',char(stor_state))
        if i_stor_state == 1
            mpc.gen_map{idx_gen, 'co2_storage_state'} = stor_state;
            mpc.gen_map{idx_gen, 'co2_storage_step'} = gen_stor_steps{i_stor_state, 'StepName'};
            mpc.gen_map{idx_gen, 'co2_storage_step_cap'} = gen_stor_steps{i_stor_state, 'step_tons'};
            mpc.gen_map{idx_gen, 'co2_storage_cost'} = gen_stor_steps{i_stor_state, 'total_cost_2013'} * mpc.gen_aux{idx_gen, 'STORAGE_CO2'};
            mpc.gen_map{idx_gen, 'co2_storage_type'} = gen_stor_steps{i_stor_state, 'stor_type'};
        else
            mpc_subset.gen_map{:, 'co2_storage_state'} = stor_state;
            mpc_subset.gen_map{:, 'co2_storage_step'} = gen_stor_steps{i_stor_state, 'StepName'};
            mpc_subset.gen_map{:, 'co2_storage_step_cap'} = gen_stor_steps{i_stor_state, 'step_tons'};
            mpc_subset.gen_map{:, 'co2_storage_cost'} = gen_stor_steps{i_stor_state, 'total_cost_2013'} * mpc.gen_aux{idx_gen, 'STORAGE_CO2'};
            mpc_subset.gen_map{:, 'co2_storage_type'} = gen_stor_steps{i_stor_state, 'stor_type'};
            [mpc, offer] = append_e4st_gen(mpc, offer, mpc_subset);
        end
    end
end

idx_ccs = ismember(mpc.gentype, {'coalccs_ret', 'ngccccs_ret', 'ngccccs_new', 'coalccs_new'});
% could be replaced by applyPTC function call, but I think inline is quicker.
% This is applied to all ccs whether new or not, because PTC resets in setup_gen
mpc.gen_aux{idx_ccs, 'PTC'} = mpc.gen_aux{idx_ccs, 'PTC'} - mpc.gen_map{idx_ccs, 'co2_storage_cost'};
mpc = updateGenCost(mpc, idx_ccs, eopt);

%% Apply Capacity Constraints:
% Retrofits
% In years after the first simyr, in case when a ccs retrofit is built,
% the original unretrofitted plant either has smaller capacity or retires, we don't want to constrain the cap of existing
% retrofits, as they already have capacity limit due to existing pmax
% So, skip existing retrofit units here, and only set constraints on new ones.
if strcmp(eopt.ccs_retro, 'T')
    idx_retro = ismember(mpc.gentype, {'coalccs_ret', 'ngccccs_ret'}) & mpc.newgen == 1;
    ids = unique(mpc.gen_map{idx_retro, 'gen_id'});
    fprintf('Applying Retrofit CCS Capacity Constraints (%d for %d total generators)\n', length(ids), sum(idx_retro))
    for id = ids'
        idx_id = mpc.gen_map{:, 'gen_id'} == id;
        
        if any(idx_id & idx_retro)
            iRetro = find(idx_id & idx_retro);
            iUnRetro = find(idx_id & ~idx_retro);
            avgcap = mpc.gen_desc{iRetro(1), 'AVG_CAP'};
            if any(idx_id & ~idx_retro)
                hr = mpc.gen_aux{iUnRetro(1), 'HR'};
            else
                hr = mpc.gen_aux{iRetro(1), 'HR'};
            end
            cappen = (ccs_retrofit_curves('CapPen', avgcap, hr) / 100);
            cap = mpc.gen(iUnRetro(1), PMAX);
            % cap = sum(mpc.gen(idx_id & idx_retro, PMAX)) / (1 - cappen); % cap max = sum of retrofits (w/o cap penalty)...
            map = double(idx_id');
            map(iRetro) = 1 / (1 - cappen);
            capMax = cap;
            capMin = -Inf;
            name = {['RetrofitCCS_', num2str(id)]};
            mpc = addCapLim(mpc, map, capMin, capMax, name);
        end
    end
end

% Decided to not limit new CCS units, since capacity was arbitrary, and raised the buildable capacity in esc.
% To eliminate any barrier to cheapest storage options being used.

% % New
% idx_ccs_new = ismember(mpc.gentype, {'coalccs_new', 'ngccccs_new'});
% ids = unique(mpc.gen_map{idx_ccs_new, 'gen_id'});
% fprintf('Applying New CCS Capacity Constraints (%d for %d total generators)\n', length(ids), sum(idx_ccs_new))
% for id = ids'
%     idx_id = mpc.gen_map{:, 'gen_id'} == id;
% 
%     map = idx_id';
%     capMax = unique(mpc.gen(idx_id, PMAX));
%     capMin = -Inf;
%     name = {['NewCCS_', num2str(id)]};
%     mpc = addCapLim(mpc, map, capMin, capMax, name);
% 
% end

%% Apply Storage Quantity Constraints
% Quantity constraints on CO2 storage options now done along with DAC in applyStorCap

%% Update MPC Dimensions
mpc = updateMPCDim(mpc);

end
