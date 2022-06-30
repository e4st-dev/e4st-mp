function [mpc, offer, eopt] = setup_gen(mpc, esc, eopt)
% setup_gen prepare set of generators for simulation

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

vfprintf(eopt.verbose, '\n\n-----Running E4ST Set-Up Module: Supply (Generators)-----\n\n')

%% Relax generator ramping constraint
% Set RAMP_10 to Inf (for all generators) to ensure that the physical ramp rate 
% never limits the change in dispatch of a unit from one hour type to another;
% 18 = RAMP_10 idx 
define_constants;
mpc.gen(:, 17:20) = 0; 
mpc.gen(:, RAMP_10) = Inf;

%% Initialize/Reset Data
mpc.gen(~ismember(mpc.genfuel, {'dl', 'dac'}), PMIN) = 0; % Initalize PMIN = 0
idx_dl = strcmp(mpc.genfuel, 'dl');
mpc.gencost = [mpc.gencost(:, 1:5), zeros(mpc.ng, 1)];
mpc.gencost(idx_dl, 1:6) = repmat([2, 0, 0, 2, 5000, 0], sum(idx_dl), 1);
mpc.gen(idx_dl, 17:20) = Inf; % Not sure if necessary...

if strcmp(eopt.update_gens, 'T')
    mpc.gen_aux{:, {'CAP_COST', 'TRANS_COST'}} = 0; % Clear annualized Capital Costs
end
mpc.gen_aux{:, {'PTC', 'ITC'}} = 0; % Clear policy subsidies
mpc.total_output = []; % Total Output Constraint
mpc.caplim = []; % Capacity Limit Constraint

% Short term Storage
mpc.short_term_storage = [(1:size(mpc.gen,1))', mpc.gen_aux{:, {'STORAGE_EFFICIENCY', 'STORAGE_E_CAPACITY'}}];

%% Initialize gen_aux columns
mpc.gen_aux{:, {'DEATHS_SO2_US_HI', 'DEATHS_NOX_US_HI', 'DEATHS_SO2_US_LO', 'DEATHS_NOX_US_LO', ...
    'DEATHS_INF_NOX', 'DEATHS_ADLT_NOX', 'DEATHS_INF_SO2', 'DEATHS_ADLT_SO2', ...
    'VSL_INF', 'VSL_ADLT', 'CH4_FUEL_CONTENT', 'DAM_CO2', 'DAM_CH4', 'CAP_CREDIT', ...
    'EMIS_CO2e', 'HIST_CF', 'TGT_AF', 'TGT_CF' 
    }} = 0;

%% Already offline
idx_offline = mpc.gen(:, GEN_STATUS) == 0;
mpc.gen(idx_offline, PG) = 0;
mpc.gen(idx_offline, PMAX) = 0;
mpc.gen(idx_offline, PMIN) = 0;

%% Remove buildable units from previous simulation that were not built to reduce size
idx_unbuilt = idx_offline & mpc.newgen == 1; % newgen not updated yet
if any(idx_unbuilt)
    mpc = removeGen(mpc, idx_unbuilt);
end

%% Update age (and vintage) of units
if strcmp(eopt.update_gens, 'T')
    mpc = setGenAge(mpc, esc.year_delta, eopt);
end

%% Initialize and set offer
offer = setGenOffer(mpc, eopt);

%% Build new generators, if applicable
if isfield(esc, 'gen_build') && strcmp(eopt.update_gens, 'T')
    [mpc, offer] = buildGen(mpc, esc, offer, eopt);
    if isfield(esc, 'custom_map')
        [mpc, esc] = addAreas(mpc, esc);
    end
end

%% Parameterize CCS: build retrofits; apply carbon transport and storage steps
if strcmp(eopt.update_gens, 'T') && (esc.year >= esc.first_build_yr) && strcmp(eopt.ccs, 'T')
    [mpc, offer] = paramCCS(mpc, esc, offer, eopt);
end

%% Parameterize DAC: apply carbon transport and storage steps
if strcmp(eopt.update_gens, 'T') && (esc.year >= esc.first_build_yr) && strcmp(eopt.dac, 'T')
    [mpc, offer] = paramDAC(mpc, esc, offer, eopt);
end

%% Set carbon storage step capacities
if strcmp(eopt.update_gens, 'T') && (esc.year >= esc.first_build_yr) && (strcmp(eopt.dac, 'T') || strcmp(eopt.ccs, 'T'))
   mpc = applyStorCap(mpc);  
end

%% Build Hydrogen Retrofits
% Co-firing retrofits
if strcmp(eopt.update_gens, 'T') && (esc.year >= esc.first_build_yr) && strcmp(eopt.hydrogen_retro, 'T')
    [mpc, offer] = paramHydrogen(mpc, esc, offer, eopt);
end

%% Set Proper Parameters for Storage
[mpc, offer] = paramStorage(mpc, esc, offer, eopt);

%% Update Generator CO2e Values
[mpc, eopt] = updateCO2e(mpc, esc, eopt);

%% Target Availability Factor / Capacity Factor for Existing Units
% Calculate historical capacity factor
mpc.gen_desc{:, 'HIST_CF'} = mpc.gen_desc{:, 'HIST_GEN'} ./ (mpc.gen_desc{:, 'HIST_CAP'} * 8760);
idx_nan = isnan(mpc.gen_desc{:, 'HIST_CF'});
mpc.gen_desc{idx_nan, 'HIST_CF'} = 0;
mpc.gen_desc{:, 'HIST_CF'} = min(mpc.gen_desc{:, 'HIST_CF'}, 1); % Cannot be greater than 1

% Non-Dispatchable: Set AFs for pre-existing units (Generators will generate up to AFs, unless curtailment)
% Wind, Solar, Hydro Run-Of-River (HYRR also set AF)
mpc.gen_aux{:, 'TGT_AF'} = 0;
idx = (mpc.newgen == min(mpc.newgen)) & (mpc.gen_desc{:, 'HIST_CF'} > 0) & (strcmp(mpc.genfuel, 'wind') | strcmp(mpc.genfuel, 'solar') | strcmp(mpc.gentype, 'hyrr'));
mpc.gen_aux{idx, 'TGT_AF'} = mpc.gen_desc{idx, 'HIST_CF'};
idx = (mpc.newgen == min(mpc.newgen)) & (mpc.gen_desc{:, 'HIST_CF'} > 0) & (strcmp(mpc.gentype, 'hyc'));
mpc.gen_aux{idx, 'TGT_AF'} = (1 + mpc.gen_desc{idx, 'HIST_CF'}) / 2;

% Dispatchable: Set CF Constraints for pre-existing units
% Biomass, Hydro Conventional and Pumped Storage
mpc.gen_aux{:, 'TGT_CF'} = 0;
idx = (mpc.newgen == min(mpc.newgen)) & (mpc.gen_desc{:, 'HIST_CF'} > 0) & (strcmp(mpc.genfuel, 'biomass') | strcmp(mpc.gentype, 'hyc') | strcmp(mpc.gentype, 'hyps'));
mpc.gen_aux{idx, 'TGT_CF'} = mpc.gen_desc{idx, 'HIST_CF'};

% Hydro Run-Of-River Constant Output set towards end of function

%% Exogenously Retire Generators (If Any)
if strcmp(eopt.update_gens, 'T')
    mpc = retireGen(mpc, offer, esc, [], eopt);
end

%% Apply capacity matching
if isfield(esc, 'cap_match')
    [mpc, offer] = matchCap(mpc, esc, offer, eopt);
end

%% Set generator availability factors
if strcmp(eopt.update_gens, 'T')
    mpc = setGenAF(mpc, esc, eopt);
    idx_wind = strcmp(mpc.genfuel(:, 1), 'wind');
    temp = sum(std(mpc.availability_factor(idx_wind,:),0,2) < 0.000000000001);
    vfprintf(eopt.verbose, '%s existing wind farms missing wind data\n', num2str(temp))
  
end

%% Target generator availability factors (ADD DL HERE TOO)
if isfield(esc, 'gen_tgt_af') && strcmp(eopt.update_gens, 'T')
    mpc = tgtGenAF(mpc, esc, eopt);
end

%% Set generator pmins
if isfield(esc, 'gen_pmin') && strcmp(eopt.update_gens, 'T')
    [mpc, offer] = setGenPmin(mpc, esc, offer, eopt);
end

%% Set generator capacity factors
if isfield(esc, 'gen_cf')
    mpc = setGenCF(mpc, esc, eopt);  
end

%% Update generator fuel costs
if isfield(esc, 'fuel_cost') && strcmp(eopt.update_gens, 'T')
    mpc = setFuelCost(mpc, esc, eopt);
end

%% Apply generator updates over time
if isfield(esc, 'gen_ann_adj') && strcmp(eopt.update_gens, 'T')
    [mpc, offer] = applyAnnAdj(mpc, offer, esc, 'gen_aux', eopt);
end

%% TEMPORARY: Transfer gencost to FOM
% https://www.eia.gov/electricity/monthly/epm_table_grapher.php?t=epmt_6_07_b
% Nuclear
idx_nuclear = strcmp(mpc.genfuel, 'nuclear');
tmp_cf = .923;
mpc.gen_aux{idx_nuclear, 'FOM'} = mpc.gen_aux{idx_nuclear, 'FOM'} + ...
    ((mpc.gen_aux{idx_nuclear, 'FUEL_COST'} .* mpc.gen_aux{idx_nuclear, 'HR'}) + ...
    mpc.gen_aux{idx_nuclear, 'VOM'} * tmp_cf);
mpc.gen_aux{idx_nuclear, {'VOM', 'FUEL_COST'}} = 0;
mpc = updateGenCost(mpc, idx_nuclear, eopt);
offer = updateOfferPrc(mpc, offer, idx_nuclear, eopt);

% Hydro Run-Of-River
idx_hyrr = strcmp(mpc.gentype, 'hyrr');
tmp_cf = 0.3718; % SNL
mpc.gen_aux{idx_hyrr, 'FOM'} = mpc.gen_aux{idx_hyrr, 'FOM'} + ...
    ((mpc.gen_aux{idx_hyrr, 'FUEL_COST'} .* mpc.gen_aux{idx_hyrr, 'HR'}) + ...
    mpc.gen_aux{idx_hyrr, 'VOM'} * tmp_cf);
mpc.gen_aux{idx_hyrr, {'VOM', 'FUEL_COST'}} = 0;
mpc = updateGenCost(mpc, idx_hyrr, eopt);
offer = updateOfferPrc(mpc, offer, idx_hyrr, eopt);

%% Past Capex
idx_new = mpc.newgen == 1;
mpc.gen_aux{idx_new, 'PAST_CAPEX'} = mpc.gen_aux{idx_new, 'CAP_COST'} + mpc.gen_aux{idx_new, 'TRANS_COST'}; % -  mpc.gen_aux{idx_new, 'ITC'};

%% Capacity Credits
cap_credit = ones(mpc.ng, 1) * .9;
cap_credit(strcmp(mpc.genfuel, 'dl')) = 0;
cap_credit(mpc.gen(:, PMIN) < 0) = 0; % DL & Charging Part of Storage
cap_credit(strcmp(mpc.genfuel, 'solar')) = .38;
cap_credit(strcmp(mpc.genfuel, 'wind')) = .13;
if any(strcmp(esc.bus_map.Properties.VariableNames, 'cap_cred_solar'))
    cap_cred_solar = getMap(mpc, esc, 'gen', 'cap_cred_solar');
    cap_credit(strcmp(mpc.genfuel, 'solar')) = cap_cred_solar{strcmp(mpc.genfuel, 'solar'), 'cap_cred_solar'};
end
if any(strcmp(esc.bus_map.Properties.VariableNames, 'cap_cred_wind'))
    cap_cred_wind = getMap(mpc, esc, 'gen', 'cap_cred_wind');
    cap_credit(strcmp(mpc.genfuel, 'wind')) = cap_cred_wind{strcmp(mpc.genfuel, 'wind'), 'cap_cred_wind'};
end
mpc.gen_aux{:, 'CAP_CREDIT'} = cap_credit;

%% Hydro Run-Of-River Constant Output
idx = (mpc.newgen == min(mpc.newgen)) & (strcmp(mpc.gentype, 'hyrr')) & (mpc.gen_desc{:, 'HIST_CF'} > 0);
mpc.availability_factor(idx, :) = repmat(mpc.gen_desc{idx, 'HIST_CF'}, 1, size(mpc.availability_factor, 2));

%% Assertions
% check generator index column in mpc.short_term_storage
check_storage = all(strcmp(mpc.genfuel(mpc.short_term_storage(:,1)',:), 'storage'));
assert(check_storage, 'Incorrect generator indices in mpc.short_term_storage')
