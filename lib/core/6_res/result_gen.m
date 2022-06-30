function res_gen = result_gen(res, mpc, offer, contab, esc, eopt)
%result_gen Summarize the results by generator

%   E4ST
%   Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).

define_constants;

%% Initialize Annual Table
res_gen.annual = table(mpc.gen(:, 1), 'VariableNames', {'bus'});

%% Locational information
if isfield(esc, 'bus_map')
    loc_map = getMap(mpc, esc, 'gen', []);
end

%% Hourly Results
hours = esc.hrs_map{:, 'hours'};
probability = esc.hrs_map{:, 'probability'};
nh = length(hours);

res_bus = res.user_result.bus_res;

res_gen.lmp_hourly = join(res_gen.annual, [array2table(mpc.bus(:, 1), 'VariableNames', {'bus'}), table(res_bus.lmp_hourly)]);
res_gen.lmp_hourly = res_gen.lmp_hourly{:, 2:end};
gen_load_hourly = join(res_gen.annual, [array2table(mpc.bus(:, 1), 'VariableNames', {'bus'}), table(res_bus.load_hourly)]);
gen_load_hourly = gen_load_hourly{:, 2:end};

res_gen.gen_hourly = zeros(size(res.base.gen, 1), nh);
res_gen.grosscs_hourly = zeros(size(res.base.gen, 1), nh);
res_gen.eleccost_hourly = zeros(size(res.base.gen, 1), nh);
res_gen.netcs_hourly = zeros(size(res.base.gen, 1), nh);
for i = 1:nh
    if i == 1
        gen = res.base.gen;
        gencost = mpc.gencost;
    else
        gen = res.cont(i - 1).gen;
        mpc_cont = apply_changes(i - 1, mpc, contab);
        gencost = mpc_cont.gencost;
    end
    res_gen.gen_hourly(:, i) = gen(:, PG);

    res_gen.grosscs_hourly(:, i) = -totcost(gencost, -gen_load_hourly(:, i));
    res_gen.eleccost_hourly(:, i) = (res_gen.lmp_hourly(:, i) .* gen_load_hourly(:, i));
    res_gen.netcs_hourly(:, i) = res_gen.grosscs_hourly(:, i) - res_gen.eleccost_hourly(:, i);
    %res_gen.totalcost_hourly(:,i) = (-totcost(gencost, -gen_load_hourly(:,i)) - (res_gen.lmp_hourly(:,i) .* gen_load_hourly(:,i))) * hours(i);

    idx_rm = ~strcmp(mpc.genfuel, 'dl'); % filter out non-dl
    res_gen.grosscs_hourly(idx_rm, i) = 0;
    res_gen.eleccost_hourly(idx_rm, i) = 0;
    res_gen.netcs_hourly(idx_rm, i) = 0;
    res_gen.totalcost_hourly(idx_rm, i) = 0;
end

%% Annual Results Table
% Descriptive and historical info
res_gen.annual{:, 'genfuel'} = mpc.genfuel;
res_gen.annual{:, 'gentype'} = mpc.gentype;
res_gen.annual{:, 'newgen'} = mpc.newgen;
res_gen.annual{:, 'status'} = mpc.gen(:, GEN_STATUS);
res_gen.annual{:, 'age'} = mpc.gen_desc{:, 'AGE'};
res_gen.annual{:, 'online_year'} = mpc.gen_desc{:, 'ON_YR'};
res_gen.annual{:, 'online_simyr'} = mpc.gen_desc{:, 'ON_SIMYR'};
res_gen.annual{:, 'retirement_year'} = mpc.gen_desc{:, 'RET_YR'};
res_gen.annual{:, 'retirement_simyr'} = mpc.gen_desc{:, 'RET_SIMYR'};
res_gen.annual{:, 'planned_retirement_year'} = mpc.gen_desc{:, 'P_RET_YR'};
res_gen.annual{:, 'retirement_reason'} = mpc.gen_desc{:, 'RET_REASON'};
res_gen.annual{:, 'invest_type'} = mpc.gen_desc{:, 'INVEST_TYPE'};
res_gen.annual{:, 'esh_type'} = mpc.gen_aux{:, 'ESH_TYPE'};
res_gen.annual{:, 'chp'} = mpc.gen_desc{:, 'CHP'};
res_gen.annual{:, 'chp_co2_multi'} = mpc.gen_aux{:, 'CHP_CO2_MULTI'};
res_gen.annual{:, 'name'} = mpc.gen_desc{:, 'NAME'};
res_gen.annual{:, 'plantid'} = mpc.gen_desc{:, 'PLANTID'};
res_gen.annual{:, 'unitid'} = mpc.gen_desc{:, 'GENID'};
res_gen.annual{:, 'hist_cap'} = mpc.gen_desc{:, 'HIST_CAP'};
res_gen.annual{:, 'hist_gen'} = mpc.gen_desc{:, 'HIST_GEN'};
res_gen.annual{:, 'hist_cf'} = res_gen.annual{:, 'hist_gen'} ./ (res_gen.annual{:, 'hist_cap'} * 8760);
res_gen.annual{:, 'hist_co2'} = res_gen.annual{:, 'hist_gen'} .* mpc.gen_aux{:, 'EMIS_CO2'};
res_gen.annual{:, 'hist_so2'} = res_gen.annual{:, 'hist_gen'} .* mpc.gen_aux{:, 'EMIS_SO2'};
res_gen.annual{:, 'hist_nox'} = res_gen.annual{:, 'hist_gen'} .* mpc.gen_aux{:, 'EMIS_NOX'};

idx_dac = strcmp(res_gen.annual{:, 'genfuel'}, 'dac');
idx_storage = strcmp(res_gen.annual{:, 'genfuel'}, 'storage');

%% Capacity
res_gen.annual{~idx_dac, 'capacity_start_mw'} = offer(~idx_dac, 2);
res_gen.annual{idx_dac, 'capacity_start_mw'} = offer(idx_dac, 4);
res_gen.annual{:, 'capacity_active_mw'} = res.reserve.qty.Rp_pos;
idx_negpmin = mpc.gen(:, PMIN) < 0; % & mpc.gen(:, PMAX) <= 0;
res_gen.annual{idx_negpmin, 'capacity_active_mw'} = res.reserve.qty.Rp_neg(idx_negpmin);
res_gen.annual{:, 'capacity_retired_mw'} = res_gen.annual{:, 'capacity_start_mw'} - res_gen.annual{:, 'capacity_active_mw'};
res_gen.annual{:, 'capacity_retired_existing_mw'} = res_gen.annual{:, 'capacity_retired_mw'} .* (res_gen.annual{:, 'newgen'} ~= 1);
res_gen.annual{:, 'capacity_retired_existing_perc'} = res_gen.annual{:, 'capacity_retired_existing_mw'} ./ res_gen.annual{:, 'capacity_start_mw'};
idx_exog_ret = strcmp(mpc.gen_desc{:, 'RET_REASON'}, 'exogenous');
res_gen.annual{:, 'capacity_retired_exog_mw'} = 0;
res_gen.annual{idx_exog_ret, 'capacity_retired_exog_mw'} = res_gen.annual{idx_exog_ret, 'capacity_retired_existing_mw'};
%idx_endog_ret = strcmp(mpc.gen_desc{:, 'RET_REASON'}, 'endogenous');
res_gen.annual{:, 'capacity_retired_endog_mw'} = 0;
res_gen.annual{~idx_exog_ret, 'capacity_retired_endog_mw'} = res_gen.annual{~idx_exog_ret, 'capacity_retired_existing_mw'};
res_gen.annual{:, 'capacity_offer_mw'} = res_gen.annual{:, 'capacity_start_mw'} - res_gen.annual{:, 'capacity_retired_exog_mw'};
res_gen.annual{:, 'capacity_invested_mw'} = res_gen.annual{:, 'capacity_active_mw'} .* (res_gen.annual{:, 'newgen'} == 1);
res_gen.annual{:, 'capacity_invested_endog_mw'} = res_gen.annual{:, 'capacity_invested_mw'} .* strcmp(mpc.gen_desc{:, 'INVEST_TYPE'}, 'standard');
res_gen.annual{:, 'capacity_invested_exog_mw'} = res_gen.annual{:, 'capacity_invested_mw'} .* (strcmp(mpc.gen_desc{:, 'INVEST_TYPE'}, 'built') | strcmp(mpc.gen_desc{:, 'INVEST_TYPE'}, 'planned'));
res_gen.annual{:, 'capacity_pmin_mw'} = mpc.gen(:, PMIN);

%% Generation
res_gen.annual{:, 'generation_mwh'} = res_gen.gen_hourly * hours;
res_gen.annual = join(res_gen.annual, res_bus.annual(:, {'bus', 'load'}));
res_gen.annual.Properties.VariableNames{end} = 'bus_load_mwh';

idx_pos = res_gen.gen_hourly > 0;
gen_hourly_pos = res_gen.gen_hourly; gen_hourly_neg = res_gen.gen_hourly;
gen_hourly_pos(~idx_pos) = 0; gen_hourly_neg(idx_pos) = 0;
res_gen.annual{:, 'generation_pos_mwh'} = gen_hourly_pos * hours;
res_gen.annual{:, 'generation_neg_mwh'} = gen_hourly_neg * hours;


res_gen.annual{:, 'capacity_factor'} = res_gen.annual{:, 'generation_mwh'} ./ (res_gen.annual{:, 'capacity_active_mw'} * sum(hours));
res_gen.annual{idx_storage, 'capacity_factor'} = res_gen.annual{idx_storage, 'generation_pos_mwh'} ./ (res_gen.annual{idx_storage, 'capacity_active_mw'} * sum(hours));
res_gen.annual{:, 'target_af'} = mpc.gen_aux{:, 'TGT_AF'};
res_gen.annual{:, 'target_cf'} = mpc.gen_aux{:, 'TGT_CF'};
res_gen.annual{:, 'min_cf_limit'} = offer(:, 13);
res_gen.annual{:, 'max_cf_limit'} = mpc.availability_factor * hours / sum(hours);
res_gen.annual{:, 'min_utilization'} = min(res_gen.gen_hourly./res_gen.annual{:, 'capacity_active_mw'}, [], 2);
res_gen.annual{:, 'max_utilization'} = max(res_gen.gen_hourly./res_gen.annual{:, 'capacity_active_mw'}, [], 2);



% idx_storage = strcmp(res_gen.annual{:, 'genfuel'}, 'storage');
% storage_gen = res_gen.annual(idx_storage,:);
% storage_gen.efficiency = storage_gen{:, 'generation_pos_mwh'}./storage_gen{:, 'generation_neg_mwh'};

%% Emissions Rates
res_gen.annual{:, 'heatrate_mmbtupermwh'} = mpc.gen_aux{:, 'HR'};
res_gen.annual{:, 'ch4_tonspermmbtu'} = mpc.gen_aux{:, 'CH4_FUEL_CONTENT'};
res_gen.annual{:, 'emisrate_co2_tonpermwh'} = mpc.gen_aux{:, 'EMIS_CO2'};
res_gen.annual{:, 'emisrate_co2e_tonpermwh'} = mpc.gen_aux{:, 'EMIS_CO2e'};
mpc2 = updateCO2e(mpc, esc, eopt, 'ch4_gwp', eopt.CH4_GWP_20y); %calculate emission rate with 20 year GWP of methane
res_gen.annual{:, 'emisrate_co2e_20y_tonpermwh'} = mpc2.gen_aux{:, 'EMIS_CO2e'};   
res_gen.annual{:, 'emisrate_ch4_tonpermwh'} = res_gen.annual{:, 'ch4_tonspermmbtu'} .* res_gen.annual{:, 'heatrate_mmbtupermwh'};
res_gen.annual{:, 'emisrate_nox_lbpermwh'} = mpc.gen_aux{:, 'EMIS_NOX'};
res_gen.annual{:, 'emisrate_so2_lbpermwh'} = mpc.gen_aux{:, 'EMIS_SO2'};
res_gen.annual{:, 'emisrate_pm25_lbpermwh'} = mpc.gen_aux{:, 'EMIS_PM25'};
res_gen.annual{:, 'storage_co2_tonpermwh'} = mpc.gen_aux{:, 'STORAGE_CO2'};

%% Emissions Totals
res_gen.annual{:, 'emis_co2_stons'} = res_gen.annual{:, 'generation_mwh'} .* res_gen.annual{:, 'emisrate_co2_tonpermwh'};
res_gen.annual{idx_storage, 'emis_co2_stons'} = res_gen.annual{idx_storage, 'generation_pos_mwh'} .* res_gen.annual{idx_storage, 'emisrate_co2_tonpermwh'};
res_gen.annual{:, 'emis_co2e_stons'} = res_gen.annual{:, 'generation_mwh'} .* res_gen.annual{:, 'emisrate_co2e_tonpermwh'};
res_gen.annual{idx_storage, 'emis_co2e_stons'} = res_gen.annual{idx_storage, 'generation_pos_mwh'} .* res_gen.annual{idx_storage, 'emisrate_co2e_tonpermwh'};
res_gen.annual{:, 'emis_co2e_20y_stons'} = res_gen.annual{:, 'generation_mwh'} .* res_gen.annual{:, 'emisrate_co2e_20y_tonpermwh'};
res_gen.annual{idx_storage, 'emis_co2e_20y_stons'} = res_gen.annual{idx_storage, 'generation_pos_mwh'} .* res_gen.annual{idx_storage, 'emisrate_co2e_20y_tonpermwh'};
res_gen.annual{:, 'emis_ch4_stons'} = res_gen.annual{:, 'generation_mwh'} .* res_gen.annual{:, 'emisrate_ch4_tonpermwh'};
res_gen.annual{idx_storage, 'emis_ch4_stons'} = res_gen.annual{idx_storage, 'generation_pos_mwh'} .* res_gen.annual{idx_storage, 'emisrate_ch4_tonpermwh'};
res_gen.annual{:, 'emis_nox_lbs'} = res_gen.annual{:, 'generation_mwh'} .* res_gen.annual{:, 'emisrate_nox_lbpermwh'};
res_gen.annual{idx_storage, 'emis_nox_lbs'} = res_gen.annual{idx_storage, 'generation_pos_mwh'} .* res_gen.annual{idx_storage, 'emisrate_nox_lbpermwh'};
res_gen.annual{:, 'emis_so2_lbs'} = res_gen.annual{:, 'generation_mwh'} .* res_gen.annual{:, 'emisrate_so2_lbpermwh'};
res_gen.annual{idx_storage, 'emis_so2_lbs'} = res_gen.annual{idx_storage, 'generation_pos_mwh'} .* res_gen.annual{idx_storage, 'emisrate_so2_lbpermwh'};
res_gen.annual{:, 'emis_pm25_lbs'} = res_gen.annual{:, 'generation_mwh'} .* res_gen.annual{:, 'emisrate_pm25_lbpermwh'};
res_gen.annual{idx_storage, 'emis_pm25_lbs'} = res_gen.annual{idx_storage, 'generation_pos_mwh'} .* res_gen.annual{idx_storage, 'emisrate_pm25_lbpermwh'};
res_gen.annual{:, 'storage_co2_stons'} = res_gen.annual{:, 'generation_mwh'} .* res_gen.annual{:, 'storage_co2_tonpermwh'};
res_gen.annual{idx_storage, 'storage_co2_stons'} = res_gen.annual{idx_storage, 'generation_pos_mwh'} .* res_gen.annual{idx_storage, 'storage_co2_tonpermwh'};

res_gen.annual{:, 'emis_co2_adj_stons'} = res_gen.annual{:, 'emis_co2_stons'} .* res_gen.annual{:, 'chp_co2_multi'};
idx_bio = strcmp(res_gen.annual{:, 'genfuel'}, 'biomass');
res_gen.annual{idx_bio, 'emis_co2_adj_stons'} = res_gen.annual{idx_bio, 'emis_co2_adj_stons'} * eopt.bio_pctCO2e; %account for the fact that not all biomass emissions are new to the atmosphere

res_gen.annual{:, 'hist_co2_adj'} = res_gen.annual{:, 'hist_co2'} .* res_gen.annual{:, 'chp_co2_multi'};
idx_bio = strcmp(res_gen.annual{:, 'genfuel'}, 'biomass');
res_gen.annual{idx_bio, 'hist_co2_adj'} = res_gen.annual{idx_bio, 'hist_co2'} * eopt.bio_pctCO2e; %account for the fact that not all biomass emissions are new to the atmosphere
%res_gen.annual{idx_bio, 'hist_co2_adj'} = 0;

%% Hourly Revenue & Costs
res_gen.annual{:, 'lmp_gen_mwh'} = sum((res_gen.gen_hourly.*hours').*res_gen.lmp_hourly, 2) ./ res_gen.annual{:, 'generation_mwh'};
idx = res_gen.annual{:, 'generation_mwh'} == 0;
res_gen.annual{idx, 'lmp_gen_mwh'} = res_gen.lmp_hourly(idx, :) * probability; % Re-compute zero generation without weighting
res_gen.annual = join(res_gen.annual, res_bus.annual(:, {'bus', 'lmp'}));
res_gen.annual.Properties.VariableNames{end} = 'lmp_bus_mwh';
res_gen.annual = join(res_gen.annual, res_bus.annual(:, {'bus', 'lmp_prob'}));
res_gen.annual.Properties.VariableNames{end} = 'lmp_bus_prob';
res_gen.annual{:, 'vom_permwh'} = mpc.gen_aux{:, 'VOM'};
res_gen.annual{:, 'fuelcost_permmbtu'} = mpc.gen_aux{:, 'FUEL_COST'};
res_gen.annual{:, 'fuelcost_permwh'} = mpc.gen_aux{:, 'FUEL_COST'} .* res_gen.annual{:, 'heatrate_mmbtupermwh'};
res_gen.annual{:, 'prod_subsidy_permwh'} = mpc.gen_aux{:, 'PTC'};
res_gen.annual{:, 'fom_permw_perhr'} = mpc.gen_aux{:, 'FOM'};
res_gen.annual{:, 'lev_capcost_permw_perhr'} = mpc.gen_aux{:, 'CAP_COST'};
res_gen.annual{:, 'lev_transcost_permw_perhr'} = mpc.gen_aux{:, 'TRANS_COST'};
res_gen.annual{:, 'routine_capex_permw_perhr'} = mpc.gen_aux{:, 'ROUTINE_CAPEX'};
res_gen.annual{:, 'past_capex_permw_perhr'} = mpc.gen_aux{:, 'PAST_CAPEX'};
% Past capex for new generators is not past. This is counted regularly in cap cost - itc.
idx_new = res_gen.annual{:, 'newgen'} == 1;
res_gen.annual{idx_new, 'past_capex_permw_perhr'} = 0;
% Portion of annual past capex that needs to be paid for generators that
% aren't new (when iterating) will vary based on when built and economic life.
% join changes the order, must use iLeft to restore original order
[mpc.gen_desc, iLeft, ~] = outerjoin(mpc.gen_desc, eopt.year_table(:,{'year','past_capex_portion'}),  'LeftKeys', 'ON_YR', 'RightKeys', 'year', 'MergeKeys', false, 'Type', 'left');
[~, sortInds] = sort(iLeft);
mpc.gen_desc = mpc.gen_desc(sortInds, :);
mpc.gen_desc{isnan(mpc.gen_desc{:,'past_capex_portion'}), 'past_capex_portion'} = 0;
res_gen.annual{:, 'past_capex_mult'} = mpc.gen_desc{:, 'past_capex_portion'};
% then multiply
res_gen.annual{:, 'past_capex_permw_perhr'} = res_gen.annual{:, 'past_capex_permw_perhr'} .* res_gen.annual{:, 'past_capex_mult'};

res_gen.annual{:, 'invest_subsidy_permw_perhr'} = mpc.gen_aux{:, 'ITC'};
res_gen.annual{:, 'fixedcost_permw_perhr'} = sum(res_gen.annual{:, {'fom_permw_perhr', 'lev_capcost_permw_perhr', 'lev_transcost_permw_perhr', 'routine_capex_permw_perhr'}}, 2);
res_gen.annual{:, 'netfixedcost_permw_perhr'} = res_gen.annual{:, 'fixedcost_permw_perhr'} - res_gen.annual{:, 'invest_subsidy_permw_perhr'};
res_gen.annual{~idx_dac, 'offerprc_permw_perhr'} = offer(~idx_dac, 1);
res_gen.annual{idx_dac, 'offerprc_permw_perhr'} = offer(idx_dac, 3);
res_gen.annual{:, 'variablecost_permwh'} = sum(res_gen.annual{:, {'vom_permwh', 'fuelcost_permwh'}}, 2);
res_gen.annual{:, 'netvariablecost_permwh'} = res_gen.annual{:, 'variablecost_permwh'} - res_gen.annual{:, 'prod_subsidy_permwh'};
res_gen.annual{:, 'gencost_permwh'} = mpc.gencost(:, 5);
res_gen.annual{:, 'prodcost_permwh'} = res_gen.annual{:, 'variablecost_permwh'} + (res_gen.annual{:, 'fixedcost_permw_perhr'} ./ res_gen.annual{:, 'capacity_factor'});
res_gen.annual{:, 'netprodcost_permwh'} = res_gen.annual{:, 'netvariablecost_permwh'} + (res_gen.annual{:, 'netfixedcost_permw_perhr'} ./ res_gen.annual{:, 'capacity_factor'});
res_gen.annual{:, 'totalcost_permwh'} = res_gen.annual{:, 'netprodcost_permwh'} + (res_gen.annual{:, 'past_capex_permw_perhr'} ./ res_gen.annual{:, 'capacity_factor'});

%% Total Revenue & Costs
% _wel variables are phased out for now.
idx_dl = strcmp(res_gen.annual{:, 'genfuel'}, 'dl');
res_gen.annual{:, 'gencost'} = res_gen.annual{:, 'gencost_permwh'} .* res_gen.annual{:, 'generation_mwh'};
res_gen.annual{idx_storage, 'gencost'} = res_gen.annual{idx_storage, 'gencost_permwh'} .* res_gen.annual{idx_storage, 'generation_pos_mwh'};
res_gen.annual{:, 'offerprc'} = res_gen.annual{:, 'offerprc_permw_perhr'} .* res_gen.annual{:, 'capacity_active_mw'} * sum(hours);
res_gen.annual{:, 'supplycost_obj'} = sum(res_gen.annual{:, {'gencost', 'offerprc'}}, 2);
res_gen.annual{idx_dl, 'supplycost_obj'} = 0;
res_gen.annual{:, 'revenue'} = res_gen.annual{:, 'lmp_gen_mwh'} .* res_gen.annual{:, 'generation_mwh'};
res_gen.annual{:, 'vom'} = res_gen.annual{:, 'vom_permwh'} .* res_gen.annual{:, 'generation_mwh'};
res_gen.annual{idx_storage, 'vom'} = res_gen.annual{idx_storage, 'vom_permwh'} .* res_gen.annual{idx_storage, 'generation_pos_mwh'};
res_gen.annual{:, 'fuelcost'} = res_gen.annual{:, 'fuelcost_permwh'} .* res_gen.annual{:, 'generation_mwh'};
res_gen.annual{idx_storage, 'fuelcost'} = res_gen.annual{idx_storage, 'fuelcost_permwh'} .* res_gen.annual{idx_storage, 'generation_pos_mwh'};
res_gen.annual{:, 'prod_subsidy'} = res_gen.annual{:, 'prod_subsidy_permwh'} .* res_gen.annual{:, 'generation_mwh'};
res_gen.annual{idx_storage, 'prod_subsidy'} = res_gen.annual{idx_storage, 'prod_subsidy_permwh'} .* res_gen.annual{idx_storage, 'generation_pos_mwh'};
res_gen.annual{:, 'fom'} = res_gen.annual{:, 'fom_permw_perhr'} .* res_gen.annual{:, 'capacity_active_mw'} * sum(hours);
res_gen.annual{:, 'capcost_obj'} = res_gen.annual{:, 'lev_capcost_permw_perhr'} .* res_gen.annual{:, 'capacity_active_mw'} * sum(hours);
res_gen.annual{:, 'capcost_obj_new'} = res_gen.annual{:, 'capcost_obj'} .* (res_gen.annual{:, 'newgen'} == 1);
%res_gen.annual{:, 'capcost_wel'} = res_gen.annual{:, 'capcost_obj'} / esc.crf / max(1, esc.year_delta);
res_gen.annual{:, 'transcost_obj'} = res_gen.annual{:, 'lev_transcost_permw_perhr'} .* res_gen.annual{:, 'capacity_active_mw'} * sum(hours);
%res_gen.annual{:, 'transcost_wel'} = res_gen.annual{:, 'transcost_obj'} / esc.crf / max(1, esc.year_delta);
res_gen.annual{:, 'routine_capex'} = res_gen.annual{:, 'routine_capex_permw_perhr'} .* res_gen.annual{:, 'capacity_active_mw'} * sum(hours);
res_gen.annual{:, 'past_capex_obj'} = res_gen.annual{:, 'past_capex_permw_perhr'} .* res_gen.annual{:, 'capacity_active_mw'} * sum(hours); % / 0.0562077; esc.ccr for 60 years
%res_gen.annual{:, 'past_capex_wel'} = res_gen.annual{:, 'past_capex_obj'} / 0.0562077 / max(1, esc.year_delta); % esc.ccr for 60 years
res_gen.annual{:, 'invest_subsidy_obj'} = res_gen.annual{:, 'invest_subsidy_permw_perhr'} .* res_gen.annual{:, 'capacity_active_mw'} * sum(hours);
%res_gen.annual{:, 'invest_subsidy_wel'} = res_gen.annual{:, 'invest_subsidy_obj'} / esc.crf / max(1, esc.year_delta);
res_gen.annual{:, 'variablecost'} = res_gen.annual{:, 'variablecost_permwh'} .* res_gen.annual{:, 'generation_mwh'};
res_gen.annual{idx_storage, 'variablecost'} = res_gen.annual{idx_storage, 'variablecost_permwh'} .* res_gen.annual{idx_storage, 'generation_pos_mwh'};
res_gen.annual{:, 'fixedcost'} = sum(res_gen.annual{:, {'fom', 'capcost_obj', 'transcost_obj', 'routine_capex'}}, 2);
res_gen.annual{:, 'prodcost'} = sum(res_gen.annual{:, {'variablecost', 'fixedcost'}}, 2);
res_gen.annual{:, 'goingforwardcost'} = sum(res_gen.annual{:, {'variablecost', 'fixedcost'}}, 2);
res_gen.annual{:, 'netvariablecost'} = res_gen.annual{:, 'variablecost'} - res_gen.annual{:, 'prod_subsidy'};
res_gen.annual{:, 'netfixedcost'} = res_gen.annual{:, 'fixedcost'} - res_gen.annual{:, 'invest_subsidy_obj'};
res_gen.annual{:, 'netprodcost'} = sum(res_gen.annual{:, {'netvariablecost', 'netfixedcost'}}, 2);


res_gen.annual{:, 20:112} = fillmissing(res_gen.annual{:, 20:112}, 'constant', 0);
%% Transfers
res_gen.annual{:, 'trans2ps'} = 0; % Payments to Producers
res_gen.annual{:, 'trans2cs'} = 0; % Payments to Consumers
res_gen.annual{:, 'trans2gov'} = 0; % Payments to Government
res_gen.annual{:, 'trans2dac'} = 0; % Payments to Direct Air Capture
res_gen.annual{:, 'trans2na'} = 0; % Payments to no one
res_gen.annual{:, 'uplift'} = 0; % Payments to Consumers
res_gen.annual{:, 'trans2seq'} = 0; % Welfare for Carbon sequesters (and transporters)

%% CO2 Storage Results
% storage cost is already incorporated into ptc, additional price paid to
% storers calced here. Don't enter if no co2 storage units.
if (strcmp(eopt.ccs, 'T') || strcmp(eopt.dac, 'T')) && any(~strcmp(mpc.gen_map{:, 'co2_storage_type'}, 'na'))
    % Cost that is used in ptc
    res_gen.annual{:, 'co2_storage_cost_total'} = res_gen.annual{:,'generation_mwh'} .* mpc.gen_map{:,'co2_storage_cost'};
    % CO2 Adjustment for leakage in EOR steps
    idx_eor = strcmp(mpc.gen_map{:, 'co2_storage_type'}, 'eor');
    res_gen.annual{idx_eor, 'emis_co2_adj_stons'} = res_gen.annual{idx_eor, 'emis_co2_adj_stons'} + eopt.eor_leakage*res_gen.annual{idx_eor, 'storage_co2_stons'};
    
    % Assuming that price paid for CO2 in each state is the price of the
    % most expensive step used in that state. 
    % Any profit earned from lower cost storage steps added to trans2seq
    
    % Find which steps are the highest used in each state
    toc.name = mpc.total_output.name;
    toc.max = res.total_output.max;
    toc.qty = res.total_output.qty;
    toc.mu = res.total_output.mu;
    toc = struct2table(toc);
    idx_step = contains(toc.name, 'CO2_Storage');
    toc = toc(idx_step, :);
    tmp = split(toc.name, '_');
    toc = [toc array2table(tmp(:, [3,4]), 'VariableNames', {'co2_storage_state', 'co2_storage_step'})];
    toc{:, 'StepNum'} = str2double(extractAfter(toc{:, 'co2_storage_step'}, 'STEP'));

    max_steps = [];
    states = unique(toc.co2_storage_state);
    for i = 1:length(states)
        state = states{i};
        % steps in each state
        idx_state = strcmp(toc.co2_storage_state, state);
        tmp = toc(idx_state, :);
        % considering less than 10 as unused.
        idx_used = tmp{:, 'qty'} > 10;
        if ~any(idx_used)
            [~, min_step] = min(tmp{:, 'StepNum'});
            max_steps = [max_steps; cell2table({state, tmp{min_step, 'co2_storage_step'}, tmp{min_step, 'StepNum'}, 0}, 'VariableNames',{'StorState' 'co2_storage_step' 'StepNum', 'shadPrc'})];  
        else
            tmp = tmp(idx_used, :);
            [~, max_step] = max(tmp{:, 'StepNum'});
            max_steps = [max_steps; cell2table({state, tmp{max_step, 'co2_storage_step'}, tmp{max_step, 'StepNum'}, tmp{max_step, 'mu'}}, 'VariableNames',{'StorState' 'co2_storage_step' 'StepNum', 'shadPrc'})];  
        end
    end
    % stor_steps
    stor_steps = esc.stor_steps;
    % get the storage state in from gen map
    res_gen.annual{:, 'ProdState'} = mpc.gen_map{:, 'state'};
    res_gen.annual{:, 'StorState'} = mpc.gen_map{:,'co2_storage_state'};
    %need to get the transport costs joined on by source state and storage state
    [~,ia] = unique(stor_steps(:, {'ProdState', 'StorState'}));
    trans_pairs = stor_steps(ia, {'ProdState', 'StorState', 'trans_cost_2013'});
    % join changes the order, must use iLeft to restore original order
    [res_gen.annual, iLeft, ~] = outerjoin(res_gen.annual, trans_pairs, 'Keys', {'ProdState', 'StorState'}, 'MergeKeys', true, 'Type', 'left');
    [~, sortInds] = sort(iLeft);
    res_gen.annual = res_gen.annual(sortInds, :);
    % then the storage prices joined onto max_steps. Right now, all steps
    % of the same number have the same storage cost
    [~, ia] = unique(stor_steps(:, 'StepNum'));
    step_costs = stor_steps(ia, {'StepNum', 'stor_cost_2013'});
    max_steps = innerjoin(max_steps, step_costs, 'Keys', 'StepNum', 'RightVariables', 'stor_cost_2013');
    % Adding on shadow prices to steps' costs
    max_steps{:, 'stor_cost_2013'} = max_steps{:, 'stor_cost_2013'} + max_steps{:, 'shadPrc'};
    % and from there onto gen results by state
    [res_gen.annual, iLeft, ~] = outerjoin(res_gen.annual, max_steps(:, {'StorState', 'stor_cost_2013'}), 'Keys', 'StorState', 'MergeKeys', true, 'Type', 'left');
    [~, sortInds] = sort(iLeft);
    res_gen.annual = res_gen.annual(sortInds, :);
    % and zero out missing values, those without storage
    res_gen.annual{ismissing(res_gen.annual{:, 'trans_cost_2013'}), 'trans_cost_2013'} = 0;
    res_gen.annual{ismissing(res_gen.annual{:, 'stor_cost_2013'}), 'stor_cost_2013'} = 0;
    % multiply to get total costs.
    res_gen.annual{:, 'co2_trans_cost'} = res_gen.annual{:, 'trans_cost_2013'} .* res_gen.annual{:,'storage_co2_stons'};
    res_gen.annual{:, 'co2_stor_price'} = res_gen.annual{:, 'stor_cost_2013'} .* res_gen.annual{:,'storage_co2_stons'};
    
    res_gen.annual{:, 'co2_storage_paid_total'} = res_gen.annual{:, 'co2_trans_cost'} + res_gen.annual{:, 'co2_stor_price'};
    
    %Drop columns used to join, they will be added back in map
    res_gen.annual = removevars(res_gen.annual, {'ProdState', 'StorState', 'trans_cost_2013', 'stor_cost_2013'});
    
    % Add to welfare transfer(so can get included in producer profit)
    res_gen.annual{:, 'trans2ps'} = res_gen.annual{:, 'trans2ps'} - res_gen.annual{:, 'co2_storage_paid_total'}; 
    res_gen.annual{:, 'trans2seq'} = res_gen.annual{:, 'co2_storage_paid_total'} - res_gen.annual{:, 'co2_storage_cost_total'};
else
    res_gen.annual{:, 'co2_storage_cost_total'} = 0;
    res_gen.annual{:, 'co2_trans_cost'} = 0;
    res_gen.annual{:, 'co2_stor_price'} = 0;
    res_gen.annual{:, 'co2_storage_paid_total'} = 0;
end

%% Policy Results
[res_gen.annual, res_gen.policy, mpc] = result_policy(res, res_gen.annual, mpc, esc, eopt);
% Adjustment of policy results for DAC welfare here.
res_gen.annual{idx_dac, 'trans2dac'} = res_gen.annual{idx_dac, 'trans2ps'} + res_gen.annual{idx_dac, 'trans2cs'};
res_gen.annual{idx_dac, 'trans2ps'} = 0;
res_gen.annual{idx_dac, 'trans2cs'} = 0;

% Adjustment of CES/RPS cost for electricity storage units. Is put into trans2cs but
% they are grouped in with the producer side. This is assuming that there
% is no other reason for storage to have anything in trans2cs.
res_gen.annual{idx_storage, 'trans2ps'} = res_gen.annual{idx_storage, 'trans2ps'} + res_gen.annual{idx_storage, 'trans2cs'};
res_gen.annual{idx_storage, 'trans2cs'} = 0;

%% Net Revenue and Transfer of Regulated Profits
%res_gen.annual{:, 'goingforwardcost'} = res_gen.annual{:, 'prodcost'} - res_gen.annual{:, 'trans2ps'};
res_gen.annual{~idx_dac, 'goingforwardcost'} = res_gen.annual{~idx_dac, 'prodcost'} - res_gen.annual{~idx_dac, 'trans2ps'};
res_gen.annual{idx_dac, 'goingforwardcost'} = res_gen.annual{idx_dac, 'prodcost'} - res_gen.annual{idx_dac, 'trans2dac'};
res_gen.annual{:, 'totalcost'} = sum(res_gen.annual{:, {'goingforwardcost', 'past_capex_obj'}}, 2);
res_gen.annual{:, 'netvarrev'} = res_gen.annual{:, 'revenue'} - res_gen.annual{:, 'netvariablecost'};
res_gen.annual{:, 'netgfrev'} = res_gen.annual{:, 'revenue'} - res_gen.annual{:, 'goingforwardcost'};
res_gen.annual{:, 'nettotalrev'} = res_gen.annual{:, 'revenue'} - res_gen.annual{:, 'totalcost'};

% Cost of Service Profits: Rebate back to consumers
% Now using reg_factor variables to determine how much of profit to rebate.
% Pre-existing units get rebated based on prelim factor. New get based on reg factor
res_gen.annual{:, 'nettotalrev_cos'} = 0;
idx_pre = res_gen.annual{:, 'newgen'} == min(res_gen.annual{:, 'newgen'}) & ~idx_dac;
idx_nonpre = res_gen.annual{:, 'newgen'} > min(res_gen.annual{:, 'newgen'}) & ~idx_dac;
res_gen.annual{idx_pre, 'nettotalrev_cos'} = res_gen.annual{idx_pre, 'nettotalrev'} .* loc_map{idx_pre, 'reg_factor_prelim'};
res_gen.annual{idx_nonpre, 'nettotalrev_cos'} = res_gen.annual{idx_nonpre, 'nettotalrev'} .* loc_map{idx_nonpre, 'reg_factor'};
res_gen.annual{:, 'trans2ps'} = res_gen.annual{:, 'trans2ps'} - res_gen.annual{:, 'nettotalrev_cos'};
res_gen.annual{:, 'trans2cs'} = res_gen.annual{:, 'trans2cs'} + res_gen.annual{:, 'nettotalrev_cos'};

res_gen.annual{:, 'goingforwardcost'} = res_gen.annual{:, 'goingforwardcost'} + res_gen.annual{:, 'nettotalrev_cos'};
res_gen.annual{:, 'totalcost'} = res_gen.annual{:, 'totalcost'} + res_gen.annual{:, 'nettotalrev_cos'};
res_gen.annual{:, 'netvarrev'} = res_gen.annual{:, 'netvarrev'} - res_gen.annual{:, 'nettotalrev_cos'};
res_gen.annual{:, 'netgfrev'} = res_gen.annual{:, 'netgfrev'} - res_gen.annual{:, 'nettotalrev_cos'};
res_gen.annual{:, 'nettotalrev'} = res_gen.annual{:, 'nettotalrev'} - res_gen.annual{:, 'nettotalrev_cos'};

%% Fuel Consumption
res_gen.annual{:, 'fuelconsumption_mmbtu'} = res_gen.annual{:, 'generation_mwh'} .* res_gen.annual{:, 'heatrate_mmbtupermwh'};
res_gen.annual{idx_storage, 'fuelconsumption_mmbtu'} = res_gen.annual{idx_storage, 'generation_pos_mwh'} .* res_gen.annual{idx_storage, 'heatrate_mmbtupermwh'};
%% Approximate Health Impacts: Deaths & Damages
% Deaths & Damages
res_gen.annual{:, 'deaths_nox_inf_perlb'} = mpc.gen_aux{:, 'DEATHS_INF_NOX'} / 2000;
res_gen.annual{:, 'deaths_nox_adlt_perlb'} = mpc.gen_aux{:, 'DEATHS_ADLT_NOX'} / 2000;
res_gen.annual{:, 'deaths_so2_inf_perlb'} = mpc.gen_aux{:, 'DEATHS_INF_SO2'} / 2000;
res_gen.annual{:, 'deaths_so2_adlt_perlb'} = mpc.gen_aux{:, 'DEATHS_ADLT_SO2'} / 2000;
res_gen.annual{:, 'damages_co2_perston'} = mpc.gen_aux{:, 'DAM_CO2'};
res_gen.annual{:, 'damages_ch4_perston'} = mpc.gen_aux{:, 'DAM_CH4'};
res_gen.annual{:, 'damages_inf_perdeath'} = mpc.gen_aux{:, 'VSL_INF'};
res_gen.annual{:, 'damages_adlt_perdeath'} = mpc.gen_aux{:, 'VSL_ADLT'};

% Deaths per mwh
res_gen.annual{:, 'deaths_nox_inf_permwh'} = res_gen.annual{:, 'deaths_nox_inf_perlb'} .* res_gen.annual{:, 'emisrate_nox_lbpermwh'};
res_gen.annual{:, 'deaths_nox_adlt_permwh'} = res_gen.annual{:, 'deaths_nox_adlt_perlb'} .* res_gen.annual{:, 'emisrate_nox_lbpermwh'};
res_gen.annual{:, 'deaths_so2_inf_permwh'} = res_gen.annual{:, 'deaths_so2_inf_perlb'} .* res_gen.annual{:, 'emisrate_so2_lbpermwh'};
res_gen.annual{:, 'deaths_so2_adlt_permwh'} = res_gen.annual{:, 'deaths_so2_adlt_perlb'} .* res_gen.annual{:, 'emisrate_so2_lbpermwh'};
res_gen.annual{:, 'deaths_inf_permwh'} = res_gen.annual{:, 'deaths_nox_inf_permwh'} + res_gen.annual{:, 'deaths_so2_inf_permwh'};
res_gen.annual{:, 'deaths_adlt_permwh'} = res_gen.annual{:, 'deaths_nox_adlt_permwh'} + res_gen.annual{:, 'deaths_so2_adlt_permwh'};
% Damages per mwh
res_gen.annual{:, 'damages_co2_permwh'} = res_gen.annual{:, 'damages_co2_perston'} .* res_gen.annual{:, 'emisrate_co2_tonpermwh'};
res_gen.annual{:, 'damages_ch4_permwh'} = res_gen.annual{:, 'damages_ch4_perston'} .* res_gen.annual{:, 'emisrate_ch4_tonpermwh'};
res_gen.annual{:, 'damages_inf_permwh'} = res_gen.annual{:, 'damages_inf_perdeath'} .* (res_gen.annual{:, 'deaths_nox_inf_permwh'} + res_gen.annual{:, 'deaths_so2_inf_permwh'});
res_gen.annual{:, 'damages_adlt_permwh'} = res_gen.annual{:, 'damages_adlt_perdeath'} .* (res_gen.annual{:, 'deaths_nox_adlt_permwh'} + res_gen.annual{:, 'deaths_so2_adlt_permwh'});
res_gen.annual{:, 'damages_permwh'} = sum(res_gen.annual{:, {'damages_co2_permwh', 'damages_ch4_permwh', 'damages_inf_permwh', 'damages_adlt_permwh'}}, 2);

% Total Deaths
res_gen.annual{:, 'deaths_inf_nox'} = res_gen.annual{:, 'deaths_nox_inf_perlb'} .* res_gen.annual{:, 'emis_nox_lbs'};
res_gen.annual{:, 'deaths_adlt_nox'} = res_gen.annual{:, 'deaths_nox_adlt_perlb'} .* res_gen.annual{:, 'emis_nox_lbs'};
res_gen.annual{:, 'deaths_nox'} = res_gen.annual{:, 'deaths_inf_nox'} + res_gen.annual{:, 'deaths_adlt_nox'};
res_gen.annual{:, 'deaths_inf_so2'} = res_gen.annual{:, 'deaths_so2_inf_perlb'} .* res_gen.annual{:, 'emis_so2_lbs'};
res_gen.annual{:, 'deaths_adlt_so2'} = res_gen.annual{:, 'deaths_so2_adlt_perlb'} .* res_gen.annual{:, 'emis_so2_lbs'};
res_gen.annual{:, 'deaths_so2'} = res_gen.annual{:, 'deaths_inf_so2'} + res_gen.annual{:, 'deaths_adlt_so2'};
res_gen.annual{:, 'deaths_total'} = res_gen.annual{:, 'deaths_nox'} + res_gen.annual{:, 'deaths_so2'};
% Total Damage of all emissions
res_gen.annual{:, 'damages_co2'} = res_gen.annual{:, 'damages_co2_perston'} .* res_gen.annual{:, 'emis_co2_adj_stons'};
res_gen.annual{:, 'damages_ch4'} = res_gen.annual{:, 'damages_ch4_perston'} .* res_gen.annual{:, 'emis_ch4_stons'};
res_gen.annual{:, 'damages_inf_nox'} = res_gen.annual{:, 'deaths_inf_nox'} .* res_gen.annual{:, 'damages_inf_perdeath'};
res_gen.annual{:, 'damages_adlt_nox'} = res_gen.annual{:, 'deaths_adlt_nox'} .* res_gen.annual{:, 'damages_adlt_perdeath'};
res_gen.annual{:, 'damages_nox'} = res_gen.annual{:, 'damages_inf_nox'} + res_gen.annual{:, 'damages_adlt_nox'};
res_gen.annual{:, 'damages_inf_so2'} = res_gen.annual{:, 'deaths_inf_so2'} .* res_gen.annual{:, 'damages_inf_perdeath'};
res_gen.annual{:, 'damages_adlt_so2'} = res_gen.annual{:, 'deaths_adlt_so2'} .* res_gen.annual{:, 'damages_adlt_perdeath'};
res_gen.annual{:, 'damages_so2'} = res_gen.annual{:, 'damages_inf_so2'} + res_gen.annual{:, 'damages_adlt_so2'};
res_gen.annual{:, 'damages_climate'} = sum(res_gen.annual{:, {'damages_co2', 'damages_ch4'}}, 2);
res_gen.annual{:, 'damages_health'} = sum(res_gen.annual{:, {'damages_nox', 'damages_so2'}}, 2);
res_gen.annual{:, 'damages_total'} = sum(res_gen.annual{:, {'damages_co2', 'damages_ch4', 'damages_nox', 'damages_so2'}}, 2);

% Low Estimate
res_gen.annual{:, 'deaths_nox_approx_lo'} = mpc.gen_aux{:, 'DEATHS_NOX_US_LO'} .* res_gen.annual{:, 'emis_nox_lbs'};
res_gen.annual{:, 'deaths_so2_approx_lo'} = mpc.gen_aux{:, 'DEATHS_SO2_US_LO'} .* res_gen.annual{:, 'emis_so2_lbs'};
res_gen.annual{:, 'deaths_total_approx_lo'} = res_gen.annual{:, 'deaths_nox_approx_lo'} + res_gen.annual{:, 'deaths_so2_approx_lo'};
res_gen.annual{:, 'damages_nox_approx_lo'} = res_gen.annual{:, 'deaths_nox_approx_lo'} .* res_gen.annual{:, 'damages_adlt_perdeath'};
res_gen.annual{:, 'damages_so2_approx_lo'} = res_gen.annual{:, 'deaths_so2_approx_lo'} .* res_gen.annual{:, 'damages_adlt_perdeath'};
res_gen.annual{:, 'damages_total_approx_lo'} = res_gen.annual{:, 'damages_nox_approx_lo'} + res_gen.annual{:, 'damages_so2_approx_lo'};
res_gen.annual{:, 'damages_total_approx_lo'} = sum(res_gen.annual{:, {'damages_co2', 'damages_ch4', 'damages_nox_approx_lo', 'damages_so2_approx_lo'}}, 2);

% High Estimate
res_gen.annual{:, 'deaths_nox_approx_hi'} = mpc.gen_aux{:, 'DEATHS_NOX_US_HI'} .* res_gen.annual{:, 'emis_nox_lbs'};
res_gen.annual{:, 'deaths_so2_approx_hi'} = mpc.gen_aux{:, 'DEATHS_SO2_US_HI'} .* res_gen.annual{:, 'emis_so2_lbs'};
res_gen.annual{:, 'deaths_total_approx_hi'} = res_gen.annual{:, 'deaths_nox_approx_hi'} + res_gen.annual{:, 'deaths_so2_approx_hi'};
res_gen.annual{:, 'damages_nox_approx_hi'} = res_gen.annual{:, 'deaths_nox_approx_hi'} .* res_gen.annual{:, 'damages_adlt_perdeath'};
res_gen.annual{:, 'damages_so2_approx_hi'} = res_gen.annual{:, 'deaths_so2_approx_hi'} .* res_gen.annual{:, 'damages_adlt_perdeath'};
res_gen.annual{:, 'damages_total_approx_hi'} = sum(res_gen.annual{:, {'damages_co2', 'damages_ch4', 'damages_nox_approx_hi', 'damages_so2_approx_hi'}}, 2);

%% Welfare (Benefit - Cost) Analysis
res_gen.annual{:, 'producerprofit_obj'} = res_gen.annual{:, 'nettotalrev'};
% This replaces these 4 components of producer profit with the _wel
% versions
% res_gen.annual{:, 'producerprofit_wel'} = res_gen.annual{:, 'producerprofit_obj'} + ...
%     res_gen.annual{:, 'capcost_obj'} - res_gen.annual{:, 'capcost_wel'} + ...
%     res_gen.annual{:, 'transcost_obj'} - res_gen.annual{:, 'transcost_wel'} + ...
%     res_gen.annual{:, 'past_capex_obj'} - res_gen.annual{:, 'past_capex_wel'} + ...
%     res_gen.annual{:, 'invest_subsidy_wel'} - res_gen.annual{:, 'invest_subsidy_obj'};

res_gen.annual{:, 'gross_cs'} = res_gen.grosscs_hourly * hours;
res_gen.annual{:, 'elec_cost'} = res_gen.eleccost_hourly * hours;
res_gen.annual{:, 'net_cs'} = res_gen.netcs_hourly * hours;
res_gen.annual{:, 'cons_payments'} = res_gen.annual{:, 'elec_cost'} - res_gen.annual{:, 'trans2cs'};

res_gen.annual{:, 'obj_val_td'} = -res.opf_results.f * sum(hours);
res_gen.annual{:, 'obj_val_bu'} = res_gen.annual{:, 'gross_cs'} - res_gen.annual{:, 'supplycost_obj'};
res_gen.annual{:, 'surplus_env'} = -res_gen.annual{:, 'damages_total'};
%res_gen.annual{:, 'surplus_consumer_td'} = res_gen.annual{:, 'obj_val_td'} - res_gen.annual{:, 'producerprofit_obj'} - res_gen.annual{:, 'surplus_merch'};
res_gen.annual{:, 'surplus_consumer_bu'} = res_gen.annual{:, 'net_cs'} + res_gen.annual{:, 'trans2cs'};
%res_gen.annual{:, 'surplus_merch'} = res_gen.annual{:, 'generation_mwh'} .* (res_gen.annual{:, 'lmp_bus_mwh'} - res_gen.annual{:, 'lmp_gen_mwh'});
%res_gen.annual{:, 'surplus_total_td'} = sum(res_gen.annual{:, {'surplus_consumer_td', 'surplus_merch', 'surplus_env', 'trans2gov', 'producerprofit_wel'}}, 2);
%res_gen.annual{:, 'surplus_total_bu'} = sum(res_gen.annual{:, {'surplus_consumer_bu', 'surplus_merch', 'surplus_env', 'trans2gov', 'producerprofit_wel'}}, 2);

%% Consumption, and correction for distribution losses
res_gen.annual{:, 'load'} = -res_gen.annual{:, 'generation_mwh'} .* strcmp(mpc.genfuel, 'dl');
res_gen.annual{:, 'load_losses'} = res_gen.annual{:, 'load'} * esc.dist_loss;
res_gen.annual{:, 'retail_sales'} = res_gen.annual{:, 'load'} * (1 - esc.dist_loss);
res_gen.annual{:, 'peak'} = max(-res_gen.gen_hourly, [], 2) .* strcmp(mpc.genfuel, 'dl');
% retail price calculated in result_byarea

if strcmp(eopt.dac, 'T')
    res_gen.annual{idx_dac, 'dac_usage_retail'} = -res_gen.annual{idx_dac, 'generation_mwh'} * (1 - esc.dist_loss);
    res_gen.annual{idx_dac, 'dac_cap_actual'} = res_gen.annual{idx_dac, 'capacity_active_mw'} * (1 - esc.dist_loss);
end

%% Locational information
if isfield(esc, 'bus_map')
    loc_map = getMap(mpc, esc, 'gen', []);
    loc_map(:, 'bus') = [];
    res_gen.annual = [res_gen.annual, loc_map];
    res_gen.annual(:, 'cn_name') = [];
    res_gen.annual{:, 'cn_name'} = mpc.gen_map{:, 'cn_name'};
    res_gen.annual(:, 'gen_id') = [];
    res_gen.annual{:, 'gen_id'} = mpc.gen_map{:, 'gen_id'};
end


%% Save Availability Factors
res_gen.afs = mpc.availability_factor;

%% Marginal Units Approximation
%idx_skip = strcmp(res_gen.annual{:,'genfuel'},'biomass') | strcmp(res_gen.annual{:,'genfuel'},'water');
%idx_marg_units = ~idx_skip & (res_gen.gen_hourly > 1) & (res_gen.gen_hourly < (res_gen.annual{:,'capacity_active_mw'} .* res_gen.afs * .99));
%res_gen.marg_units = ones(size(idx_marg_units)) .* idx_marg_units;
