function user_res = result_byarea(user_res, esc, areas, eopt)
% result_byarea aggregates E4ST results into area level

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Areas
%all_areas = {'grid','interconnection','nation','state','nerc','subnerc','iso'};
all_areas = {'grid', 'interconnection', 'nation', 'state'};
%all_areas = {'state'};
for area = areas
    if ~any(strcmp(area, all_areas))
        all_areas = [all_areas, area];
    end
end
rpt_areas = areas; % Save for end
areas = all_areas;

%%
% Rep Hour Names
hours = esc.hrs_map{:, 'hours'};
probability = esc.hrs_map{:, 'probability'};
repHrs = {};
for i = 1:length(hours)
    repHrs{end + 1} = strcat('C', num2str(i));
end  

%% Bus
if isfield(user_res, 'bus_res')
    ann_res = user_res.bus_res.annual;
    vars_sum = {'load', 'peak', 'load_losses', 'retail_sales', 'uncurtailed_load', 'curtailed_load'};
    vars_wgt = {'lmp', 'lmp_prob', 'latitude', 'longitude'};
    

    count = 0;
    for area = areas
        vars_group = area;
        %cur_table = unique(ann_res(:, vars_group));

        if ~any(strcmp(area, ann_res.Properties.VariableNames))
            continue
        end

        count = count + 1;

        tmp = varfun(@sum, ann_res, 'GroupingVariables', vars_group, 'InputVariables', vars_sum);
        tmp.Properties.VariableNames(length(vars_group)+2:end) = vars_sum;
        sum_table = [array2table(repmat(area, height(tmp), 1)), tmp];
        sum_table.Properties.VariableNames(1:2) = {'area', 'subarea'};

        tmp = ann_res(:, [vars_group, vars_wgt]);
        tmp{:, vars_wgt} = tmp{:, vars_wgt} .* ann_res{:, 'load'};
        tmp = varfun(@sum, tmp, 'GroupingVariables', vars_group, 'InputVariables', vars_wgt);
        tmp = tmp(:, (end -(length(vars_wgt) - 1)):end);
        tmp.Properties.VariableNames = vars_wgt;
        tmp{:, vars_wgt} = tmp{:, vars_wgt} ./ sum_table{:, 'load'};
        sum_table = [sum_table, tmp];

        if count == 1
            user_res.area_res.bus = sum_table;
        else
            user_res.area_res.bus = [user_res.area_res.bus; sum_table];
        end
    end
end

%% Generators
if isfield(user_res, 'gen_res')
    ann_res = user_res.gen_res.annual;

    vars_sum = {'hist_cap', 'hist_gen', 'hist_co2', 'hist_co2_adj', 'hist_so2', 'hist_nox', 'capacity_start_mw', 'capacity_offer_mw', 'capacity_active_mw', ...
        'capacity_retired_mw', 'capacity_retired_existing_mw', 'capacity_retired_endog_mw', 'capacity_retired_exog_mw', ...
        'capacity_invested_mw', 'capacity_invested_endog_mw', 'capacity_invested_exog_mw', ...
        'generation_mwh', 'generation_pos_mwh', 'generation_neg_mwh', 'emis_co2_stons', 'emis_co2_adj_stons', 'emis_co2e_stons', 'emis_co2e_20y_stons', 'emis_ch4_stons', 'emis_nox_lbs', 'emis_so2_lbs', 'emis_pm25_lbs', 'storage_co2_stons', ...
        'revenue', 'vom', 'fuelcost', ...
        'gencost', 'offerprc', 'supplycost_obj', 'gross_cs', 'elec_cost', 'net_cs', ...
        'fom', 'capcost_obj', 'capcost_obj_new', 'transcost_obj', ...
        'routine_capex', 'past_capex_obj', 'variablecost', 'fixedcost', 'goingforwardcost', ...
        'totalcost', 'netvarrev', 'netgfrev', 'nettotalrev', 'fuelconsumption_mmbtu', ...
        'deaths_inf_nox', 'deaths_adlt_nox', 'deaths_nox', 'deaths_inf_so2', 'deaths_adlt_so2', 'deaths_so2', 'deaths_total', ...
        'damages_co2', 'damages_ch4', 'damages_inf_nox', 'damages_adlt_nox', 'damages_nox', 'damages_inf_so2', 'damages_adlt_so2', 'damages_so2', 'damages_total', ...
        'damages_health', 'damages_climate', ...
        'producerprofit_obj',  ...
        'netprodcost', 'netfixedcost', 'netvariablecost', 'prodcost', 'prod_subsidy', 'invest_subsidy_obj',  ...
        'trans2cs', 'trans2ps', 'trans2gov', 'trans2na', 'uplift', 'cons_payments', 'nettotalrev_cos', ...
        'load', 'load_losses', 'retail_sales', 'peak', ... %'retail_prc', ...
        'deaths_so2_approx_lo', 'deaths_nox_approx_lo', 'deaths_total_approx_lo', 'deaths_so2_approx_hi', 'deaths_nox_approx_hi', 'deaths_total_approx_hi', ...
        'damages_so2_approx_lo', 'damages_nox_approx_lo', 'damages_total_approx_lo', 'damages_so2_approx_hi', 'damages_nox_approx_hi', 'damages_total_approx_hi', ...
        'obj_val_bu', 'surplus_env', 'surplus_consumer_bu', ...
        'co2_storage_cost_total', 'co2_trans_cost', 'co2_stor_price', 'co2_storage_paid_total', 'trans2seq'};
    % 'surplus_merch', , 'surplus_consumer_bu', 'surplus_consumer_td', 'surplus_total_td', 'surplus_total_bu',
    % 'producerprofit_wel', 'invest_subsidy_wel', 'past_capex_wel', 'transcost_wel', 'capcost_wel',}

    if strcmp(eopt.dac, 'T')
        vars_sum = [vars_sum, {'trans2dac', 'dac_usage_retail', 'dac_cap_actual'}];
    end

    vars_wgt = {'hist_cf', 'obj_val_td', 'max_cf_limit', 'min_cf_limit', 'max_utilization', 'min_utilization', 'capacity_factor', 'capacity_retired_existing_perc', ...
        'heatrate_mmbtupermwh', 'ch4_tonspermmbtu', 'emisrate_co2_tonpermwh', 'emisrate_co2e_tonpermwh', 'emisrate_ch4_tonpermwh', ...
        'emisrate_nox_lbpermwh', 'emisrate_so2_lbpermwh', 'storage_co2_tonpermwh', ...
        'lmp_bus_prob', 'lmp_bus_mwh', 'lmp_gen_mwh', 'vom_permwh', 'fuelcost_permmbtu', 'fuelcost_permwh', 'fom_permw_perhr', 'lev_capcost_permw_perhr', 'lev_transcost_permw_perhr', ...
        'routine_capex_permw_perhr', 'past_capex_permw_perhr', 'fixedcost_permw_perhr', ...
        'netfixedcost_permw_perhr', 'offerprc_permw_perhr', 'variablecost_permwh', 'gencost_permwh', 'totalcost_permwh', ...
        'netprodcost_permwh', 'prodcost_permwh', 'netvariablecost_permwh', 'prod_subsidy_permwh', 'invest_subsidy_permw_perhr', ...
        'deaths_nox_inf_perlb', 'deaths_nox_adlt_perlb', 'deaths_so2_inf_perlb', 'deaths_so2_adlt_perlb', 'damages_co2_perston', 'damages_ch4_perston', 'damages_inf_perdeath', 'damages_adlt_perdeath', ...
        'deaths_nox_inf_permwh', 'deaths_nox_adlt_permwh', 'deaths_so2_inf_permwh', 'deaths_so2_adlt_permwh', 'deaths_inf_permwh', 'deaths_adlt_permwh', ...
        'damages_co2_permwh', 'damages_ch4_permwh', 'damages_inf_permwh', 'damages_adlt_permwh', 'damages_permwh', ...
        'latitude', 'longitude'};

    % Keep same order
    vars_all = [];
    for var = ann_res.Properties.VariableNames
        if any(strcmp(var, vars_sum)) || any(strcmp(var, vars_wgt))
            vars_all = [vars_all, var];
        end
    end

    vars_status = {'on'}; %, 'ret'

    idx_dl = strcmp(ann_res{:, 'genfuel'}, 'dl');
    for status = vars_status
        switch char(status)
            case 'on'
                idx_gen = ann_res{:, 'status'} == 1 | ann_res{:, 'retirement_year'} == esc.year;
            case 'ret'
                idx_gen = ann_res{:, 'status'} == 0;
        end

        count = 0;
        for area = areas

            % Aggregate results for generators, by area
            vars_group = area;
            tmp_ann_res = ann_res(idx_gen & ~idx_dl, :);
            %cur_table = unique(tmp_ann_res(:, vars_group));

            if ~any(strcmp(area, tmp_ann_res.Properties.VariableNames))
                continue
            end

            count = count + 1;

            tmp = varfun(@sum, tmp_ann_res, 'GroupingVariables', vars_group, 'InputVariables', vars_sum);
            tmp.Properties.VariableNames(length(vars_group)+2:end) = vars_sum;
            sum_table = [array2table(repmat(area, height(tmp), 1)), tmp];
            sum_table.Properties.VariableNames(1:2) = {'area', 'subarea'};

            tmp = tmp_ann_res(:, [vars_group, vars_wgt]);
            tmp{:, vars_wgt} = tmp{:, vars_wgt} .* tmp_ann_res{:, 'capacity_start_mw'};
            tmp = varfun(@sum, tmp, 'GroupingVariables', vars_group, 'InputVariables', vars_wgt);
            tmp = tmp(:, (end -(length(vars_wgt) - 1)):end);
            tmp.Properties.VariableNames = vars_wgt;
            tmp{:, vars_wgt} = tmp{:, vars_wgt} ./ sum_table{:, 'capacity_start_mw'};
            sum_table = [sum_table, tmp];
            sum_table = [sum_table(:, 1:(length(vars_group) + 2)), sum_table(:, vars_all)];
            sum_table{:, 'capacity_factor'} = sum_table{:, 'generation_mwh'} ./ (sum_table{:, 'capacity_active_mw'} * sum(hours));

            if count == 1
                user_res.area_res.(strjoin([{'gen'}, status], {'_'})) = sum_table;
            else
                user_res.area_res.(strjoin([{'gen'}, status], {'_'})) = [user_res.area_res.(strjoin([{'gen'}, status], {'_'})); sum_table];
            end

            % Aggregate results for generators by area and genfuel
            tmp_ann_res = ann_res(idx_gen, :);
            vars_group = [area, {'genfuel'}];
            %cur_table = unique(tmp_ann_res(:, vars_group));

            tmp = varfun(@sum, tmp_ann_res, 'GroupingVariables', vars_group, 'InputVariables', vars_sum);
            tmp.Properties.VariableNames(length(vars_group)+2:end) = vars_sum;
            sum_table = [array2table(repmat(area, height(tmp), 1)), tmp];
            sum_table.Properties.VariableNames(1:2) = {'area', 'subarea'};

            tmp = tmp_ann_res(:, [vars_group, vars_wgt]);
            tmp{:, vars_wgt} = tmp{:, vars_wgt} .* tmp_ann_res{:, 'capacity_start_mw'};
            tmp = varfun(@sum, tmp, 'GroupingVariables', vars_group, 'InputVariables', vars_wgt);
            tmp = tmp(:, (end -(length(vars_wgt) - 1)):end);
            tmp.Properties.VariableNames = vars_wgt;
            tmp{:, vars_wgt} = tmp{:, vars_wgt} ./ sum_table{:, 'capacity_start_mw'};
            sum_table = [sum_table, tmp];
            sum_table = [sum_table(:, 1:(length(vars_group) + 2)), sum_table(:, vars_all)];
            sum_table{:, 'capacity_factor'} = sum_table{:, 'generation_mwh'} ./ (sum_table{:, 'capacity_active_mw'} * sum(hours));

            if count == 1
                user_res.area_res.(strjoin([{'gen'}, status, vars_group(2:end)], {'_'})) = sum_table;
            else
                user_res.area_res.(strjoin([{'gen'}, status, vars_group(2:end)], {'_'})) = [user_res.area_res.(strjoin([{'gen'}, status, vars_group(2:end)], {'_'})); sum_table];
            end

            % Aggregate results for generators by area and gentype
            vars_group = [area, {'genfuel'}, {'gentype'}];
            cur_table = unique(tmp_ann_res(:, vars_group));

            tmp = varfun(@sum, tmp_ann_res, 'GroupingVariables', vars_group, 'InputVariables', vars_sum);
            tmp.Properties.VariableNames(length(vars_group)+2:end) = vars_sum;
            sum_table = [array2table(repmat(area, height(tmp), 1)), tmp];
            sum_table.Properties.VariableNames(1:2) = {'area', 'subarea'};

            tmp = tmp_ann_res(:, [vars_group, vars_wgt]);
            tmp{:, vars_wgt} = tmp{:, vars_wgt} .* tmp_ann_res{:, 'capacity_start_mw'};
            tmp = varfun(@sum, tmp, 'GroupingVariables', vars_group, 'InputVariables', vars_wgt);
            tmp = tmp(:, (end -(length(vars_wgt) - 1)):end);
            tmp.Properties.VariableNames = vars_wgt;
            tmp{:, vars_wgt} = tmp{:, vars_wgt} ./ sum_table{:, 'capacity_start_mw'};
            sum_table = [sum_table, tmp];
            sum_table = [sum_table(:, 1:(length(vars_group) + 2)), sum_table(:, vars_all)];
            sum_table{:, 'capacity_factor'} = sum_table{:, 'generation_mwh'} ./ (sum_table{:, 'capacity_active_mw'} * sum(hours));

            if count == 1
                user_res.area_res.(strjoin([{'gen'}, status, {'gentype'}], {'_'})) = sum_table;
            else
                user_res.area_res.(strjoin([{'gen'}, status, {'gentype'}], {'_'})) = [user_res.area_res.(strjoin([{'gen'}, status, {'gentype'}], {'_'})); sum_table];
            end

            % Marginal Units by Gentype
            if isfield(user_res.gen_res, 'marg_units')
                idx_marg_units = user_res.gen_res.marg_units(idx_gen, :);
                if size(idx_marg_units, 2) < length(repHrs)
                    idx_marg_units = [idx_marg_units, zeros(size(idx_marg_units, 1), length(repHrs) - size(idx_marg_units, 2))];
                end

                var = 'generation_mwh';
                op = @sum;
                tmp_res = tmp_ann_res(:, vars_group);
                tmp_res{:, 'VarName'} = {var};
                tmp_res = [tmp_res, array2table(idx_marg_units .* tmp_ann_res{:, var}, 'VariableNames', repHrs)];
                tmp_hr_res = varfun(op, tmp_res, 'GroupingVariables', [vars_group, {'VarName'}], 'InputVariables', repHrs);
                tmp_hr_res.Properties.VariableNames((end-length(repHrs)+1):end) = repHrs;
                tmp_hr_res{:, 'annual'} = tmp_hr_res{:, repHrs} * hours; % Hour-weighted total
                sum_table = [array2table(repmat(area, height(tmp_hr_res), 1)), tmp_hr_res(:, [vars_group, {'VarName', 'annual'}, repHrs])];
                sum_table.Properties.VariableNames(1:2) = {'area', 'subarea'};

                vars = {'heatrate_mmbtupermwh', 'emisrate_co2_tonpermwh', 'emisrate_nox_lbpermwh', 'emisrate_so2_lbpermwh'};
                op = @mean;
                for var = vars
                    tmp_res = tmp_ann_res(:, vars_group);
                    tmp_res{:, 'VarName'} = var;
                    tmp_res = [tmp_res, array2table(idx_marg_units .* tmp_ann_res{:, var}, 'VariableNames', repHrs)];
                    tmp_hr_res = varfun(op, tmp_res, 'GroupingVariables', [vars_group, {'VarName'}], 'InputVariables', repHrs);
                    tmp_hr_res.Properties.VariableNames((end-length(repHrs)+1):end) = repHrs;
                    tmp_hr_res{:, 'annual'} = tmp_hr_res{:, repHrs} * probability; % Hour-weighted average
                    tmp_sum_table = [array2table(repmat(area, height(tmp_hr_res), 1)), tmp_hr_res(:, [vars_group, {'VarName', 'annual'}, repHrs])];
                    tmp_sum_table.Properties.VariableNames(1:2) = {'area', 'subarea'};
                    sum_table = [sum_table; tmp_sum_table];
                end

                if count == 1
                    user_res.area_res.(strjoin([{'gen'}, status, {'marg_unit'}], {'_'})) = sum_table;
            else
                    user_res.area_res.(strjoin([{'gen'}, status, {'marg_unit'}], {'_'})) = [user_res.area_res.(strjoin([{'gen'}, status, {'marg_unit'}], {'_'})); sum_table];
                end
            end

        end
    end
end

%% Hourly Generation / Prices
if isfield(user_res, 'gen_res')
    hrs = esc.hrs_map{:, 'rep_hr'}';
    gen_hrs = [];
    lmp_hrs = [];
    for hr = hrs
        gen_hrs = [gen_hrs, {strjoin([{'GEN_'}, hr], '')}];
        lmp_hrs = [lmp_hrs, {strjoin([{'LMP_'}, hr], '')}];
    end
    ann_res = user_res.gen_res.annual;
    gen_hourly = array2table(user_res.gen_res.gen_hourly, 'VariableNames', gen_hrs);
    lmp_hourly = array2table(user_res.gen_res.lmp_hourly, 'VariableNames', lmp_hrs);

    hourly = [ann_res(:, {'genfuel', 'gentype', 'capacity_active_mw', 'generation_mwh', 'capacity_factor', 'lmp_gen_mwh'}), gen_hourly, lmp_hourly];
    hourly{:, lmp_hrs} = hourly{:, lmp_hrs} .* hourly{:, gen_hrs}; % For weighted average
    hourly{:, {'capacity_factor', 'lmp_gen_mwh'}} = hourly{:, {'capacity_factor', 'lmp_gen_mwh'}} .* hourly{:, 'generation_mwh'}; % For weighted average

    vars_sum = [{'capacity_active_mw', 'generation_mwh'}, gen_hrs];
    vars_wgt = [{'capacity_factor', 'lmp_gen_mwh'}, lmp_hrs];
    vars_all = hourly.Properties.VariableNames(3:end);

    count = 0;
    for area = areas
        tmp_res = [ann_res(:, area), hourly];

        if ~any(strcmp(area, tmp_res.Properties.VariableNames))
            continue
        end

        count = count + 1;

        vars_group = tmp_res.Properties.VariableNames(1:3);

        tmp = varfun(@sum, tmp_res, 'GroupingVariables', vars_group, 'InputVariables', vars_sum);
        tmp.Properties.VariableNames(length(vars_group)+2:end) = vars_sum;
        sum_table = [array2table(repmat(area, height(tmp), 1)), tmp];
        sum_table.Properties.VariableNames(1:2) = {'area', 'subarea'};

        % LMP
        tmp = tmp_res(:, [vars_group, vars_wgt]);
        tmp = varfun(@sum, tmp, 'GroupingVariables', vars_group, 'InputVariables', vars_wgt);
        tmp = tmp(:, (end -(length(vars_wgt) - 1)):end);
        tmp.Properties.VariableNames = vars_wgt;
        tmp{:, lmp_hrs} = tmp{:, lmp_hrs} ./ sum_table{:, gen_hrs};
        tmp{:, {'capacity_factor', 'lmp_gen_mwh'}} = tmp{:, {'capacity_factor', 'lmp_gen_mwh'}} ./ sum_table{:, 'generation_mwh'};
        sum_table = [sum_table, tmp];
        sum_table = [sum_table(:, 1:(length(vars_group) + 2)), sum_table(:, vars_all)];
        sum_table{:, 'capacity_factor'} = sum_table{:, 'generation_mwh'} ./ (sum_table{:, 'capacity_active_mw'} * sum(hours));

        if count == 1
            user_res.area_res.(strjoin([{'gen'}, status, {'gentype_hourly'}], {'_'})) = sum_table;
        else
            user_res.area_res.(strjoin([{'gen'}, status, {'gentype_hourly'}], {'_'})) = [user_res.area_res.(strjoin([{'gen'}, status, {'gentype_hourly'}], {'_'})); sum_table];
        end
    end
end

%% Flows
if isfield(user_res, 'branch_res')
    vars_sum = {'flow_mwh', 'merch_surplus'};
    vars_sum_norev = {'merch_surplus'};

    flow_list = [];
    % allows for no dc lines.
    if isfield(user_res.branch_res, 'annual_dc') 
        types = {'ac', 'dc'};
    else
        types = {'ac'};
    end
    
    for type = types
        %flow_list = [];
        for area = areas
            ann_res = user_res.branch_res.(['annual_', char(type)]);

            if ~any(strcmp(area, esc.bus_map.Properties.VariableNames))
                continue
            end

            map = esc.bus_map(:, [{'bus'}, area]);
            ann_res{:, 'bus'} = ann_res{:, 'fbus'};
            ann_res = join(ann_res, map);
            ann_res.Properties.VariableNames(area) = {'farea'};
            ann_res{:, 'bus'} = ann_res{:, 'tbus'};
            ann_res = join(ann_res, map);
            ann_res.Properties.VariableNames(area) = {'tarea'};

            hrly_flow = user_res.branch_res.(['flow_', char(type), '_hourly']);
            if size(hrly_flow, 2) < length(repHrs)
                hrly_flow = [hrly_flow, zeros(size(hrly_flow, 1), length(repHrs) - size(hrly_flow, 2))];
            end
            ann_res = [ann_res, array2table(hrly_flow, 'VariableNames', repHrs(1:size(hrly_flow, 2)))];

            subareas = unique(esc.bus_map{:, area});
            %row_names = array2table(subareas, 'VariableNames', {'farea'});
            %flow_table = array2table(zeros(length(subareas), length(subareas)), 'VariableNames', subareas, 'RowNames', subareas');

            for subarea = subareas'
                idx = strcmp(ann_res{:, 'farea'}, subarea);
                tmp_ann_res = ann_res(idx, :);
                tmp = varfun(@sum, tmp_ann_res, 'GroupingVariables', {'tarea', 'br_type'}, 'InputVariables', [vars_sum, repHrs(1:size(hrly_flow, 2))]);
                tmp.Properties.VariableNames(4:end) = [vars_sum, repHrs(1:size(hrly_flow, 2))];

                idx = strcmp(ann_res{:, 'tarea'}, subarea);
                tmp_ann_res = ann_res(idx, :);
                tmp_rev = varfun(@sum, tmp_ann_res, 'GroupingVariables', {'farea', 'br_type'}, 'InputVariables', [vars_sum, repHrs(1:size(hrly_flow, 2))]);
                tmp_rev.Properties.VariableNames(4:end) = [vars_sum, repHrs(1:size(hrly_flow, 2))];

                tmp_rev{:, 4:end} = -1 * tmp_rev{:, 4:end};
                tmp_rev.Properties.VariableNames(1) = {'tarea'};
                tmp_rev{:, vars_sum_norev} = -1 * tmp_rev{:, vars_sum_norev}; % flip back sign; sum regardless of direction
                if isempty(tmp) && isempty(tmp_rev)
                    continue
                end
                tmp_final = [tmp; tmp_rev];
                %tmp_final.Properties.VariableNames(3:end) = [vars_sum repHrs(1:size(hrly_flow, 2))];
                tmp_final = varfun(@sum, tmp_final, 'GroupingVariables', {'tarea', 'br_type'}, 'InputVariables', [{'GroupCount'}, vars_sum, repHrs(1:size(hrly_flow, 2))]);
                tmp_final.Properties.VariableNames(3:end) = [{'Count_Area', 'Count_Branch'}, vars_sum, repHrs(1:size(hrly_flow, 2))];
                %tmp_rows = outerjoin(row_names, tmp, 'Type', 'Left', 'LeftKeys', {'farea'}, 'RightKeys', {'tarea'});
                %flow_table{:, subarea} = tmp_rows{:, 'flow_mwh'};
                flow_list = [flow_list; [array2table(repmat([area, subarea], height(tmp_final), 1), 'VariableNames', {'area', 'farea'}), tmp_final]];
            end

            %user_res.area_res.branch.(char(type)).(char(area)) = flow_table;
        end
        %user_res.area_res.(['branch_' (char(type))]) = flow_list;
    end
    user_res.area_res.branch = flow_list;
end

%% Update Bus Results
% New variables added here need to be added to list in result_user_agg
% This part is for transferring area results from gen, branch, and elsewhere to bus area results
% And using them to calculate new area results.
ur_bus = user_res.area_res.bus;

% initialize ur_bus columns that will be added
ur_bus{:, 'ms_int'} = 0;
ur_bus{:, 'ms_ext'} = 0;
ur_bus{:, 'ms_tot'} = 0;
ur_bus{:, 'trans2cs_tot'} = 0;
ur_bus{:, 'surplus_env'} = 0;
ur_bus{:, 'surplus_prod'} = 0;
ur_bus{:, 'trans2cs_gen'} = 0;
ur_bus{:, 'surplus_gov'} = 0;
ur_bus{:, 'elec_cost'} = 0;
ur_bus{:, 'cons_payments'} = 0;
ur_bus{:, 'dist_payments'} = 0;
ur_bus{:, 'uplift'} = 0;
ur_bus{:, 'retail_prc'} = 0;
ur_bus{:, 'surplus_cons'} = 0;
ur_bus{:, 'surplus_total'} = 0;
for area = areas
    idx_bus_area = strcmp(area, ur_bus{:, 'area'});
    subareas = unique(ur_bus{idx_bus_area, 'subarea'})';
    for subarea = subareas
        idx_bus = idx_bus_area & strcmp(subarea, ur_bus{:, 'subarea'});

        % Uplift
        ur_gen = user_res.area_res.gen_on;
        idx_gen = strcmp(ur_gen{:, 'area'}, area) & strcmp(ur_gen{:, 'subarea'}, subarea);
        if ~any(idx_gen)
            gen_vars = {'trans2cs_gen', 'surplus_env', 'surplus_prod', 'surplus_gov'};
            ur_bus{idx_bus, gen_vars} = 0;
        else
            ur_gen_fuel = user_res.area_res.gen_on_genfuel;
            idx_gen_fuel = strcmp(ur_gen_fuel{:, 'area'}, area) & strcmp(ur_gen_fuel{:, 'subarea'}, subarea);
            ur_bus{idx_bus, 'trans2cs_gen'} = sum(ur_gen_fuel{idx_gen_fuel, 'trans2cs'});
            %gen_vars = {'capacity_start_mw', 'capacity_active_mw', 'generation_mwh', ...
            %            'emis_co2_stons','emis_co2e_stons','emis_so2_lbs','emis_nox_lbs'}; %'trans2cs',
            %ur_bus{idx_bus, gen_vars} = ur_gen{idx_gens, gen_vars};

            % Welfare
            ur_bus{idx_bus, 'surplus_env'} = ur_gen{idx_gen, 'surplus_env'};
            ur_bus{idx_bus, 'surplus_prod'} = ur_gen{idx_gen, 'producerprofit_obj'};
            % Warning, this only sums trans2gov from generators. If there
            % is any trans2gov from dl, this will not pick it up! Summing
            % from gen_on_genfuel will pick it up.
            ur_bus{idx_bus, 'surplus_gov'} = ur_gen{idx_gen, 'trans2gov'};
        end

        ur_br = user_res.area_res.branch;
        idx_br = strcmp(ur_br{:, 'area'}, area) & ...
            strcmp(ur_br{:, 'farea'}, subarea) & ...
            strcmp(ur_br{:, 'tarea'}, subarea); % internal
        ms_int = sum(ur_br{idx_br, 'merch_surplus'}) / 2;
        idx_br = strcmp(ur_br{:, 'area'}, area) & ~idx_br & ...
            (strcmp(ur_br{:, 'farea'}, subarea) | ...
            strcmp(ur_br{:, 'tarea'}, subarea)); % external
        ms_ext = sum(ur_br{idx_br, 'merch_surplus'}) / 4;
        ms_tot = ms_int + ms_ext;
        ur_bus{idx_bus, 'ms_int'} = ms_int;
        ur_bus{idx_bus, 'ms_ext'} = ms_ext;
        ur_bus{idx_bus, 'ms_tot'} = ms_tot;

        ur_bus{idx_bus, 'trans2cs_tot'} = ur_bus{idx_bus, 'trans2cs_gen'} + ms_tot;
        
        ur_bus{idx_bus, 'elec_cost'} = ur_bus{idx_bus, 'load'} * ur_bus{idx_bus, 'lmp'};
        ur_bus{idx_bus, 'cons_payments'} = ur_bus{idx_bus, 'elec_cost'} - ur_bus{idx_bus, 'trans2cs_tot'};
        ur_bus{idx_bus, 'dist_payments'} = ur_bus{idx_bus, 'retail_sales'} * esc.dist_cost;
        ur_bus{idx_bus, 'uplift'} = -ur_bus{idx_bus, 'trans2cs_tot'} / ur_bus{idx_bus, 'load'};
        ur_bus{idx_bus, 'retail_prc'} = sum(ur_bus{idx_bus, {'lmp', 'uplift'}}) / (1 - esc.dist_loss) + esc.dist_cost;

        % Consumer Surplus
        e = esc.elasticity;
        p1 = 17811;
        p0 = ur_bus{idx_bus, 'retail_prc'};
        q0 = ur_bus{idx_bus, 'uncurtailed_load'} * (1 - esc.dist_loss); % Uncurtailed retail sales
        A = q0 / (p0^e);
        cs = (A / (e + 1)) * (p1^(e + 1) - p0^(e + 1));
        cs_tot = cs - ((ur_bus{idx_bus, 'curtailed_load'} * (1 - esc.dist_loss)) * ((5000 * (1 + esc.dist_loss)) + esc.dist_cost));
        ur_bus{idx_bus, 'surplus_cons'} = cs_tot;

        ur_bus{idx_bus, 'surplus_total'} = sum(ur_bus{idx_bus, {'surplus_env', 'surplus_prod', 'surplus_gov', 'surplus_cons'}}, 2);
    end
end
user_res.area_res.bus = ur_bus;

%% Areas
%all_areas = {'grid','interconnection','nation','state','nerc','subnerc','iso'};
area_fields = fieldnames(user_res.area_res)';
for field = area_fields
    data = user_res.area_res.(char(field));
    idx_keep = zeros(height(data), 1);
    if strcmp(field, 'bus')
        continue
    else
        for area = rpt_areas
            idx_keep = idx_keep | strcmp(data{:, 'area'}, area);
        end
        user_res.area_res.(char(field)) = data(idx_keep, :);
    end
end
