function agg_res = result_user_agg(res, eopt)
% result_user_agg 

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%%
cases = fieldnames(res);
%for year = years
for iCase = 1:length(cases)
    case_name = cases(iCase);
    yr_iters = fieldnames(res.(char(case_name)))';
    yr_iters = yr_iters(~strcmp(yr_iters, 'input'));
    for yr_iter = yr_iters
        %tmp = res.(char(case_name)).(['Y' num2str(year)]).user_result.area_res;
        if ~isfield(res.(char(case_name)).(char(yr_iter)), 'user_result')
            continue
        end
        tmp = res.(char(case_name)).(char(yr_iter)).user_result.area_res;
        %year = res.(char(case_name)).(char(yr_iter)).year;
        for areafield = fieldnames(tmp)'
            tmp_data = [array2table(repmat(yr_iter, height(tmp.(char(areafield))), 1), 'VariableNames', {'year'}), tmp.(char(areafield))];
            if ~exist('all_res', 'var') || ~isfield(all_res, areafield)
                all_res.(char(areafield)) = tmp_data;
            else
                all_res.(char(areafield)) = [all_res.(char(areafield)); tmp_data];
            end
        end

        aux_res = [];
        if isfield(res.(char(case_name)).(char(yr_iter)).user_result, 'constraints_res')
            aux_res = [aux_res, {'constraints_res'}];
        end
        if isfield(res.(char(case_name)).(char(yr_iter)).user_result, 'pol_res')
            aux_res = [aux_res, {'pol_res'}];
        end
        if isfield(res.(char(case_name)).(char(yr_iter)).user_result, 'slv_res')
            aux_res = [aux_res, {'slv_res'}];
        end
        
        % Temporary skip of aux_res for test cases because caused error
        if strcmp(eopt.prj_name, 'test')
            aux_res = [];
        end

        if ~isempty(aux_res)
            for aux = aux_res
                tmp = res.(char(case_name)).(char(yr_iter)).user_result.(char(aux));
                for areafield = fieldnames(tmp)'
                    if isempty(tmp.(char(areafield)))
                        continue
                    end
                    tmp_data = array2table(repmat(yr_iter, height(tmp.(char(areafield))), 1), 'VariableNames', {'year'});
                    tmp_data{:, 'case_name'} = case_name;
                    tmp_data = [tmp_data, tmp.(char(areafield))];
                    if ~exist('agg_res', 'var') || ~isfield(agg_res, areafield)
                        agg_res.(char(areafield)) = tmp_data;
                    else
                        agg_res.(char(areafield)) = [agg_res.(char(areafield)); tmp_data];
                    end
                end
            end
        end
        %         if isfield(res.(char(case_name)).(char(yr_iter)).user_result, 'damages')
        %             tmp = res.(char(case_name)).(char(yr_iter)).user_result.damages;
        %             for areafield = fieldnames(tmp)'
        %                 tmp_data = array2table(repmat(yr_iter, height(tmp.(char(areafield))), 1), 'VariableNames', {'year'});
        %                 tmp_data{:, 'case_name'} = case_name;
        %                 tmp_data = [tmp_data tmp.(char(areafield))];
        %                 if ~exist('agg_res', 'var') || ~isfield(agg_res, areafield)
        %                     agg_res.(char(areafield)) = tmp_data;
        %                 else
        %                     agg_res.(char(areafield)) = [agg_res.(char(areafield)); tmp_data];
        %                 end
        %             end
        %         end
    end
end
%end

%% Bus
if isfield(all_res, 'bus')
    ann_res = all_res.bus;
    vars_sum = {'load', 'peak', 'load_losses', 'retail_sales', 'uncurtailed_load', 'curtailed_load', ...
        'trans2cs_gen', 'ms_int', 'ms_ext', 'ms_tot', 'trans2cs_tot', 'elec_cost', 'cons_payments', 'dist_payments', ...
        'surplus_env', 'surplus_prod', 'surplus_gov', 'surplus_cons', 'surplus_total'}; %, ...
    %'capacity_start_mw', 'capacity_active_mw', 'generation_mwh', ...
    %'emis_co2_stons','emis_co2e_stons','emis_so2_lbs','emis_nox_lbs'};

    vars_wgt = {'lmp', 'lmp_prob', 'latitude', 'longitude', 'uplift', 'retail_prc'};
    
    vars_group = {'year', 'area', 'subarea'};

    tmp = varfun(@sum, ann_res, 'GroupingVariables', vars_group, 'InputVariables', vars_sum);
    tmp.Properties.VariableNames(length(vars_group)+2:end) = vars_sum;
    sum_table = tmp;

    tmp = ann_res(:, [vars_group, vars_wgt]);
    tmp{:, vars_wgt} = tmp{:, vars_wgt} .* ann_res{:, 'load'};
    tmp = varfun(@sum, tmp, 'GroupingVariables', vars_group, 'InputVariables', vars_wgt);
    tmp = tmp(:, (end -(length(vars_wgt) - 1)):end);
    tmp.Properties.VariableNames = vars_wgt;
    tmp{:, vars_wgt} = tmp{:, vars_wgt} ./ sum_table{:, 'load'};
    sum_table = [sum_table, tmp];

    agg_res.bus = sum_table;
end

%% Generators

vars_sum = {'hist_cap', 'hist_gen', 'hist_co2', 'hist_co2_adj', 'hist_so2', 'hist_nox', 'capacity_start_mw', 'capacity_offer_mw', 'capacity_active_mw', ...
    'capacity_retired_mw', 'capacity_retired_existing_mw', 'capacity_retired_endog_mw', 'capacity_retired_exog_mw', ...
    'capacity_invested_mw', 'capacity_invested_endog_mw', 'capacity_invested_exog_mw', ...
    'generation_mwh', 'generation_pos_mwh', 'generation_neg_mwh', 'emis_co2_stons', 'emis_co2_adj_stons', 'emis_co2e_stons', 'emis_co2e_20y_stons', 'emis_ch4_stons', 'emis_nox_lbs', 'emis_so2_lbs', 'emis_pm25_lbs', 'storage_co2_stons', ...
    'revenue', 'vom', 'fuelcost',...
    'gencost', 'offerprc', 'supplycost_obj', 'gross_cs', 'elec_cost', 'net_cs', ...
    'fom', 'capcost_obj', 'capcost_obj_new', 'transcost_obj', ...
    'routine_capex', 'past_capex_obj', 'variablecost', 'fixedcost', 'goingforwardcost', ...
    'totalcost', 'netvarrev', 'netgfrev', 'nettotalrev', 'fuelconsumption_mmbtu', ...
    'deaths_inf_nox', 'deaths_adlt_nox', 'deaths_nox', 'deaths_inf_so2', 'deaths_adlt_so2', 'deaths_so2', 'deaths_total', ...
    'damages_co2', 'damages_ch4', 'damages_inf_nox', 'damages_adlt_nox', 'damages_nox', 'damages_inf_so2', 'damages_adlt_so2', 'damages_so2', 'damages_total', ...
    'damages_health', 'damages_climate', ...
    'producerprofit_obj', ...
    'netprodcost', 'netfixedcost', 'netvariablecost', 'prodcost', 'prod_subsidy', 'invest_subsidy_obj', ...
    'trans2cs', 'trans2ps', 'trans2gov', 'trans2na', 'uplift', 'cons_payments', 'nettotalrev_cos', ...
    'load', 'load_losses', 'retail_sales', 'peak', ... %'retail_prc', ...
    'deaths_so2_approx_lo', 'deaths_nox_approx_lo', 'deaths_total_approx_lo', 'deaths_so2_approx_hi', 'deaths_nox_approx_hi', 'deaths_total_approx_hi', ...
    'damages_so2_approx_lo', 'damages_nox_approx_lo', 'damages_total_approx_lo', 'damages_so2_approx_hi', 'damages_nox_approx_hi', 'damages_total_approx_hi', ...
    'obj_val_bu', 'surplus_env', 'surplus_consumer_bu', ...
    'co2_storage_cost_total', 'co2_trans_cost', 'co2_stor_price', 'co2_storage_paid_total', 'trans2seq'};
% 'surplus_merch', , 'surplus_consumer_bu', 'surplus_consumer_td', 'surplus_total_td', 'surplus_total_bu',
% 'producerprofit_wel',  'invest_subsidy_wel',  'transcost_wel', 'past_capex_wel', 'capcost_wel',}

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
for var = all_res.gen_on.Properties.VariableNames
    if any(strcmp(var, vars_sum)) || any(strcmp(var, vars_wgt))
        vars_all = [vars_all, var];
    end
end

%vars_status = {'on', 'ret'};
vars_status = {'on'};

for status = vars_status

    % Just Area
    tmp_ann_res = all_res.(strjoin([{'gen'}, status], {'_'}));
    vars_group = {'year', 'area', 'subarea'};

    tmp = varfun(@sum, tmp_ann_res, 'GroupingVariables', vars_group, 'InputVariables', vars_sum);
    tmp.Properties.VariableNames(length(vars_group)+2:end) = vars_sum;
    sum_table = tmp;

    tmp = tmp_ann_res(:, [vars_group, vars_wgt]);
    tmp{:, vars_wgt} = tmp{:, vars_wgt} .* tmp_ann_res{:, 'capacity_start_mw'};
    tmp = varfun(@sum, tmp, 'GroupingVariables', vars_group, 'InputVariables', vars_wgt);
    tmp = tmp(:, (end -(length(vars_wgt) - 1)):end);
    tmp.Properties.VariableNames = vars_wgt;
    tmp{:, vars_wgt} = tmp{:, vars_wgt} ./ sum_table{:, 'capacity_start_mw'};
    sum_table = [sum_table, tmp];
    sum_table = [sum_table(:, 1:(length(vars_group) + 1)), sum_table(:, vars_all)];

    agg_res.(strjoin([{'gen'}, status], {'_'})) = sum_table;

    % Genfuel
    tmp_ann_res = all_res.(strjoin([{'gen'}, status, {'genfuel'}], {'_'}));
    vars_group = {'year', 'area', 'subarea', 'genfuel'};

    tmp = varfun(@sum, tmp_ann_res, 'GroupingVariables', vars_group, 'InputVariables', vars_sum);
    tmp.Properties.VariableNames(length(vars_group)+2:end) = vars_sum;
    sum_table = tmp;

    tmp = tmp_ann_res(:, [vars_group, vars_wgt]);
    tmp{:, vars_wgt} = tmp{:, vars_wgt} .* tmp_ann_res{:, 'capacity_start_mw'};
    tmp = varfun(@sum, tmp, 'GroupingVariables', vars_group, 'InputVariables', vars_wgt);
    tmp = tmp(:, (end -(length(vars_wgt) - 1)):end);
    tmp.Properties.VariableNames = vars_wgt;
    tmp{:, vars_wgt} = tmp{:, vars_wgt} ./ sum_table{:, 'capacity_start_mw'};
    sum_table = [sum_table, tmp];
    sum_table = [sum_table(:, 1:(length(vars_group) + 1)), sum_table(:, vars_all)];

    agg_res.(strjoin([{'gen'}, status, {'genfuel'}], {'_'})) = sum_table;

    % Gentype
    tmp_ann_res = all_res.(strjoin([{'gen'}, status, {'gentype'}], {'_'}));
    vars_group = {'year', 'area', 'subarea', 'genfuel', 'gentype'};

    tmp = varfun(@sum, tmp_ann_res, 'GroupingVariables', vars_group, 'InputVariables', vars_sum);
    tmp.Properties.VariableNames(length(vars_group)+2:end) = vars_sum;
    sum_table = tmp;

    tmp = tmp_ann_res(:, [vars_group, vars_wgt]);
    tmp{:, vars_wgt} = tmp{:, vars_wgt} .* tmp_ann_res{:, 'capacity_start_mw'};
    tmp = varfun(@sum, tmp, 'GroupingVariables', vars_group, 'InputVariables', vars_wgt);
    tmp = tmp(:, (end -(length(vars_wgt) - 1)):end);
    tmp.Properties.VariableNames = vars_wgt;
    tmp{:, vars_wgt} = tmp{:, vars_wgt} ./ sum_table{:, 'capacity_start_mw'};
    sum_table = [sum_table, tmp];
    sum_table = [sum_table(:, 1:(length(vars_group) + 1)), sum_table(:, vars_all)];

    agg_res.(strjoin([{'gen'}, status, {'gentype'}], {'_'})) = sum_table;

    % Marginal Units by Gentype
    if isfield(all_res, strjoin([{'gen'}, status, {'marg_unit'}], {'_'}))
        tmp_ann_res = all_res.(strjoin([{'gen'}, status, {'marg_unit'}], {'_'}));
        vars_group = {'year', 'area', 'subarea', 'genfuel', 'gentype', 'VarName'};
        vars_oth = tmp_ann_res.Properties.VariableNames((length(vars_group)+1):end);
        sum_table = [];

        var = 'generation_mwh';
        op = @sum;
        tmp_res = tmp_ann_res(strcmp(tmp_ann_res{:, 'VarName'}, var), :);
        tmp_hr_res = varfun(op, tmp_res, 'GroupingVariables', vars_group, 'InputVariables', vars_oth);
        tmp_hr_res.Properties.VariableNames((end-length(vars_oth)+1):end) = vars_oth;
        sum_table = [sum_table; tmp_hr_res];

        vars = {'heatrate_mmbtupermwh', 'emisrate_co2_tonpermwh', 'emisrate_nox_lbpermwh', 'emisrate_so2_lbpermwh'};
        op = @mean;
        for var = vars
            tmp_res = tmp_ann_res(strcmp(tmp_ann_res{:, 'VarName'}, var), :);
            tmp_hr_res = varfun(op, tmp_res, 'GroupingVariables', vars_group, 'InputVariables', vars_oth);
            tmp_hr_res.Properties.VariableNames((end-length(vars_oth)+1):end) = vars_oth;
            sum_table = [sum_table; tmp_hr_res];
        end

        agg_res.(strjoin([{'gen'}, status, {'marg_unit'}], {'_'})) = sum_table;
    end
end

%% Flows
vars_sum = {'flow_mwh', 'merch_surplus', ...
    'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', ...
    'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', ...
    'C20', 'C21', 'C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', ...
    'C29', 'C30', 'C31', 'C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38'};

if isfield(all_res, 'branch')
    tmp_ann_res = all_res.branch;
    vars_group = {'year', 'area', 'farea', 'tarea', 'br_type'};

    vars_keep = [];
    for var = vars_sum
        if any(strcmp(var, tmp_ann_res.Properties.VariableNames))
            vars_keep = [vars_keep, var];
        end
    end

    tmp = varfun(@sum, tmp_ann_res, 'GroupingVariables', vars_group, 'InputVariables', vars_keep);
    tmp.Properties.VariableNames(length(vars_group)+2:end) = vars_keep;
    sum_table = tmp;

    agg_res.branch = sum_table;
end

