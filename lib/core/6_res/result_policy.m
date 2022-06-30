function [res_gen, res_pol, mpc] = result_policy(res, res_gen, mpc, esc, eopt)
% result_policy processes results of power sector policies

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

warning('off', 'MATLAB:table:RowsAddedNewVars');

%% Initialize Table

info = esc.policy;
info = filterInfo(info, esc);

info.value = info.value(:, ['Y', num2str(esc.year)]);

if isfield(info, 'joint')
    info_list = info.joint{:, :};
else
    info_list = (1:height(info.value));
end

res_pol = [];

% Things which have to pay for credits for elec they consume (demand side)
idx_ps_dl = ismember(info.genfuel{:, :}, {'dl', 'dac', 'storage'}) & ismember(info.pol{:, 'type'}, ...
    {'rps', 'ces', 'ces_coal', 'ces_ngcc'});

% For elec consumers, this section brings the pctage they have to pay into gen_wgt
info.gen_wgt{idx_ps_dl, :} = info.value{idx_ps_dl, :};

%% Calculate Policy Results
for i = unique(info_list)'
    idx_info = i == info_list;

    cur_info = filterStruct(info, idx_info);
    [idx_gen, gen_wgt] = getInfoIdx(mpc, esc, cur_info, 'gen', 1);

    if ~any(idx_gen)
        continue
    end

    pol_type = unique(cur_info.pol{:, 'type'});
    pol_name = unique(cur_info.pol{:, 'name'});
    pol_val = unique(cur_info.value{:, :});
    f_wel = strjoin([{'trans2'}, unique(cur_info.wel{:, 'f_wel'})], '');
    t_wel = strjoin([{'trans2'}, unique(cur_info.wel{:, 't_wel'})], '');
    %strjoin([{num2str(iter_num)} iter_name], '_');

    mpc.gen_map{:, [char(pol_name), '_', char(pol_type)]} = {'out'};
    mpc.gen_map{idx_gen, [char(pol_name), '_', char(pol_type)]} = {'in'};

    switch char(pol_type)
        case 'calib_qty'
        case 'calib_prc'
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_val, (pol_qty * pol_val)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];
        case 'ptc'
            res_gen{idx_gen, [char(pol_name), '_val']} = pol_val;
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_val']} .* res_gen{idx_gen, 'generation_mwh'};
            res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
            res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

            pol_qty = sum(res_gen{idx_gen, 'generation_mwh'});
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_val, (pol_qty * pol_val)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];

        case 'itc'
            res_gen{idx_gen, [char(pol_name), '_val']} = pol_val;
            res_gen{idx_gen, [char(pol_name), '_cost_obj']} = res_gen{idx_gen, [char(pol_name), '_val']} .* res_gen{idx_gen, 'capcost_obj'};
            res_gen{idx_gen, [char(pol_name), '_cost_wel']} = res_gen{idx_gen, [char(pol_name), '_cost_obj']} / esc.crf / max(1, esc.year_delta);
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_cost_obj']};
            res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
            res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

            pol_qty = sum(res_gen{idx_gen, 'capacity_invested_mw'});
            pol_cost = sum(res_gen{idx_gen, [char(pol_name), '_cost_obj']});
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_val, pol_cost], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];

        case 'emis_prc_co2'
            res_gen{idx_gen, [char(pol_name), '_mu']} = pol_val;
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_mu']} .* res_gen{idx_gen, 'emis_co2_stons'};
            res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
            res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

            pol_qty = sum(res_gen{idx_gen, 'emis_co2_stons'});
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_val, (pol_qty * pol_val)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];

        case 'emis_prc_co2e'
            res_gen{idx_gen, [char(pol_name), '_mu']} = pol_val;
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_mu']} .* res_gen{idx_gen, 'emis_co2e_stons'};
            res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
            res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

            pol_qty = sum(res_gen{idx_gen, 'emis_co2e_stons'});
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_val, (pol_qty * pol_val)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];

        case 'emis_prc_co2e_20y'
            res_gen{idx_gen, [char(pol_name), '_mu']} = pol_val;
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_mu']} .* res_gen{idx_gen, 'emis_co2e_20y_stons'};
            res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
            res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

            pol_qty = sum(res_gen{idx_gen, 'emis_co2e_20y_stons'});
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_val, (pol_qty * pol_val)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];
        
        case 'emis_prc_nox'
            res_gen{idx_gen, [char(pol_name), '_mu']} = pol_val;
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_mu']} .* res_gen{idx_gen, 'emis_nox_lbs'};
            res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
            res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

            pol_qty = sum(res_gen{idx_gen, 'emis_nox_lbs'});
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_val, (pol_qty * pol_val)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];

        case 'emis_prc_so2'
            res_gen{idx_gen, [char(pol_name), '_mu']} = pol_val;
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_mu']} .* res_gen{idx_gen, 'emis_so2_lbs'};
            res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
            res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

            pol_qty = sum(res_gen{idx_gen, 'emis_so2_lbs'});
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_val, (pol_qty * pol_val)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];

        case 'emis_prc_pm'
            res_gen{idx_gen, [char(pol_name), '_mu']} = pol_val;
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_mu']} .* res_gen{idx_gen, 'emis_pm25_lbs'};
            res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
            res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

            pol_qty = sum(res_gen{idx_gen, 'emis_pm25_lbs'});
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_val, (pol_qty * pol_val)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];
        
        case 'emis_cap_co2'
            res_gen{idx_gen, [char(pol_name), '_cap']} = pol_val;
            idx_mu = strcmp(mpc.total_output.name, char(pol_name));
            res_gen{idx_gen, [char(pol_name), '_mu']} = res.total_output.mu(idx_mu);
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_mu']} .* res_gen{idx_gen, 'emis_co2_stons'};
            res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
            res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

            pol_prc = res.total_output.mu(idx_mu);
            pol_qty = sum(res_gen{idx_gen, 'emis_co2_stons'});
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_prc, (pol_qty * pol_prc)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];

        case 'emis_cap_co2e'
            res_gen{idx_gen, [char(pol_name), '_cap']} = pol_val;
            idx_mu = strcmp(mpc.total_output.name, char(pol_name));
            res_gen{idx_gen, [char(pol_name), '_mu']} = res.total_output.mu(idx_mu);
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_mu']} .* res_gen{idx_gen, 'emis_co2e_stons'};
            res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
            res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

            pol_prc = res.total_output.mu(idx_mu);
            pol_qty = sum(res_gen{idx_gen, 'emis_co2e_stons'});
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_prc, (pol_qty * pol_prc)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];

        case 'emis_cap_nox'
            res_gen{idx_gen, [char(pol_name), '_cap']} = pol_val;
            idx_mu = strcmp(mpc.total_output.name, char(pol_name));
            res_gen{idx_gen, [char(pol_name), '_mu']} = res.total_output.mu(idx_mu);
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_mu']} .* res_gen{idx_gen, 'emis_nox_lbs'};
            res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
            res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

            pol_prc = res.total_output.mu(idx_mu);
            pol_qty = sum(res_gen{idx_gen, 'emis_nox_lbs'});
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_prc, (pol_qty * pol_prc)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];

        case 'emis_cap_so2'
            res_gen{idx_gen, [char(pol_name), '_cap']} = pol_val;
            idx_mu = strcmp(mpc.total_output.name, char(pol_name));
            res_gen{idx_gen, [char(pol_name), '_mu']} = res.total_output.mu(idx_mu);
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_mu']} .* res_gen{idx_gen, 'emis_so2_lbs'};
            res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
            res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

            pol_prc = res.total_output.mu(idx_mu);
            pol_qty = sum(res_gen{idx_gen, 'emis_so2_lbs'});
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_prc, (pol_qty * pol_prc)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];

        case 'co2_storage_prc'
            res_gen{idx_gen, [char(pol_name), '_mu']} = pol_val;
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_mu']} .* res_gen{idx_gen, 'storage_co2_stons'};
            res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
            res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

            pol_qty = sum(res_gen{idx_gen, 'storage_co2_stons'});
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_val, (pol_qty * pol_val)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];

        case 'co2_storage_cap'
            res_gen{idx_gen, [char(pol_name), '_cap']} = pol_val;
            idx_mu = strcmp(mpc.total_output.name, char(pol_name));
            res_gen{idx_gen, [char(pol_name), '_mu']} = res.total_output.mu(idx_mu);
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_mu']} .* res_gen{idx_gen, 'storage_co2_stons'};
            res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
            res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

            pol_prc = res.total_output.mu(idx_mu);
            pol_qty = sum(res_gen{idx_gen, 'storage_co2_stons'});
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_prc, (pol_qty * pol_prc)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];

        case 'irm'
            if ~isempty(mpc.caplim)
                idx_mu = strcmp(mpc.caplim.name, char(pol_name));
                if any(idx_mu)
                    map = mpc.caplim.map(idx_mu, :)';
                    idx_gen = map > 0;

                    %mpc.gen_map{:, [char(pol_name) '_' char(pol_type)]} = ones(size(idx_gen, 1), 1) .* map;
                    mpc.gen_map{:, [char(pol_name), '_', char(pol_type)]} = {'out'};
                    mpc.gen_map{idx_gen, [char(pol_name), '_', char(pol_type)]} = {'in'};

                    res_gen{idx_gen, [char(pol_name), '_req']} = mpc.caplim.min(idx_mu);
                    res_gen{idx_gen, [char(pol_name), '_mu']} = -res.caplim.mu(idx_mu);
                    res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_mu']} .* res_gen{idx_gen, 'capacity_active_mw'} .* map(idx_gen) * sum(esc.hrs_map{:, 'hours'});
                    res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
                    res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

                    pol_qty = sum(res_gen{idx_gen, 'capacity_active_mw'}.*map(idx_gen));
                    pol_prc = -res.caplim.mu(idx_mu);
                    pol_cost = sum(res_gen{idx_gen, [char(pol_name), '_cost']});
                    res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_prc, pol_cost], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];
                end
            end

        case 'uplift'
            res_gen{idx_gen, [char(pol_name), '_mu']} = pol_val;
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_mu']} .* -res_gen{idx_gen, 'generation_mwh'};
            res_gen{idx_gen, 'uplift'} = res_gen{idx_gen, 'uplift'} + res_gen{idx_gen, [char(pol_name), '_cost']};
            %res_gen{idx_gen, 'trans2ps'} = res_gen{idx_gen, 'trans2ps'} + res_gen{idx_gen, [char(pol_name) '_cost']};
            %res_gen{idx_gen, 'trans2cs'} = res_gen{idx_gen, 'trans2cs'} - res_gen{idx_gen, [char(pol_name) '_cost']};

            pol_qty = sum(res_gen{idx_gen, 'generation_mwh'});
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_val, (pol_qty * pol_val)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];

        case 'rps'
            run resultPS;
        case 'rps_qty'
            if pol_val >= 0 && pol_val <= 1 % target value, convert to fraction
                idx_area_bus = getInfoIdx(mpc, esc, rmfield(cur_info, {'genfuel', 'gentype'}), 'bus');
                pol_val = sum(mpc.ann_load(idx_area_bus)) * pol_val;
            end
            res_gen{idx_gen, [char(pol_name), '_req']} = pol_val;
            idx_mu = strcmp(mpc.total_output.name, char(pol_name));
            res_gen{idx_gen, [char(pol_name), '_mu']} = res.total_output.mu(idx_mu);
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_mu']} .* res_gen{idx_gen, 'generation_mwh'};
            res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
            res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

            pol_qty = sum(res_gen{idx_gen, 'generation_mwh'});
            pol_prc = res.total_output.mu(idx_mu);
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_prc, (pol_qty * pol_prc)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];

        case 'rps_qty_max'
            if pol_val >= 0 && pol_val <= 1 % target value, convert to fraction
                idx_area_bus = getInfoIdx(mpc, esc, rmfield(cur_info, {'genfuel', 'gentype'}), 'bus');
                pol_val = sum(mpc.ann_load(idx_area_bus)) * pol_val;
            end
            res_gen{idx_gen, [char(pol_name), '_req']} = pol_val;
            idx_mu = strcmp(mpc.total_output.name, char(pol_name));
            res_gen{idx_gen, [char(pol_name), '_mu']} = res.total_output.mu(idx_mu);
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_mu']} .* res_gen{idx_gen, 'generation_mwh'};
            res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
            res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

            pol_qty = sum(res_gen{idx_gen, 'generation_mwh'});
            pol_prc = res.total_output.mu(idx_mu);
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_prc, (pol_qty * pol_prc)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];

        case 'rps_prc'
            res_gen{idx_gen, [char(pol_name), '_mu']} = pol_val;

            res_gen{idx_gen, [char(pol_name), '_credits']} = res_gen{idx_gen, 'generation_mwh'};
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_mu']} .* res_gen{idx_gen, [char(pol_name), '_credits']};
            res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
            res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

            pol_qty = sum(res_gen{idx_gen, [char(pol_name), '_credits']});
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_val, (pol_qty * pol_val)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];

        case 'ces'
            run resultPS;
        case 'ces_coal'
            run resultPS;
        case 'clean_futures_act_2021'
            run resultPS;            
        case 'ces_ngcc'
            run resultPS;

        case 'ces_prc'
            res_gen{idx_gen, [char(pol_name), '_mu']} = pol_val;

            gen_er = res_gen{:, 'emisrate_co2e_tonpermwh'};
            credit = 1 - (gen_er / bm_er);
            credit(strcmp(mpc.genfuel, 'biomass')) = min(.5*(2204 / 2000), credit(strcmp(mpc.genfuel, 'biomass')));
            credit(credit < 0) = 0;

            res_gen{idx_gen, [char(pol_name), '_credits']} = res_gen{idx_gen, 'generation_mwh'} .* credit(idx_gen);
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_mu']} .* res_gen{idx_gen, [char(pol_name), '_credits']};
            res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
            res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

            pol_qty = sum(res_gen{idx_gen, [char(pol_name), '_credits']});
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_val, (pol_qty * pol_val)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];

        case 'ces_prc_ngcc'
            res_gen{idx_gen, [char(pol_name), '_mu']} = pol_val;

            bm_er = .4 * (2204 / 2000);

            gen_er = res_gen{:, 'emisrate_co2e_tonpermwh'};

            credit = 1 - (gen_er / bm_er);
            credit(strcmp(mpc.genfuel, 'biomass')) = min(.5*(2204 / 2000), credit(strcmp(mpc.genfuel, 'biomass')));
            credit(strcmp(mpc.gentype, 'coal_cofire')) = .5 * .15;
            credit(credit < 0) = 0;
            credit(strcmp(mpc.genfuel, 'dac')) = - gen_er(strcmp(mpc.genfuel, 'dac')) / bm_er;

            res_gen{idx_gen, [char(pol_name), '_credits']} = res_gen{idx_gen, 'generation_mwh'} .* credit(idx_gen);
            res_gen{idx_gen, [char(pol_name), '_cost']} = res_gen{idx_gen, [char(pol_name), '_mu']} .* res_gen{idx_gen, [char(pol_name), '_credits']};
            res_gen{idx_gen, f_wel} = res_gen{idx_gen, f_wel} - res_gen{idx_gen, [char(pol_name), '_cost']};
            res_gen{idx_gen, t_wel} = res_gen{idx_gen, t_wel} + res_gen{idx_gen, [char(pol_name), '_cost']};

            pol_qty = sum(res_gen{idx_gen, [char(pol_name), '_credits']});
            res_pol = [res_pol; [table(pol_type), table(pol_name), array2table([pol_val, pol_qty, pol_val, (pol_qty * pol_val)], 'VariableNames', {'pol_val', 'qty', 'prc', 'cost'})]];


        otherwise
            %vfprintf(eopt.verbose, 'Policy type %s not recognized\n', char(pol_type))
            fprintf('Policy type %s not recognized\n', char(pol_type))
    end
end

