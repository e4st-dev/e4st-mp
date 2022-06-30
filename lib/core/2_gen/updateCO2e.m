function [mpc, eopt] = updateCO2e(mpc, esc, eopt, varargin)
%UPDATE_CO2e This function reads in important parameters from the esc.gen_ann_adj
%sheet that control CO2-equivalent emissions. These parameters include:
% methane leakage rates of coal and natural gas units
% global warming potential of methane
% the leakage rate of CO2 using enhanced oil recovery
% the fraction of biomass emissions that are new to the atmosphere
%Using these values, the function updates the generator CO2e emissions. 
% It adds important parameters to the options file for later use.
%
%INPUTS
%
%mpc (struct)... matpower case file
%esc (struct)... e4st case file
%eopt (struct)... e4st options structure
%
% OPTIONAL NAMED ARGUMENTS
%
% These can be specified in the function call as
% name-value pairs. They override the default values in esc file, but are
% not saved in eopt for later use.
%
% ch4_gwp (double)     ... the global warming potential of methane (stons CO2/stons CH4)
% bio_pctCO2e (double) ... percentage of biomass emissions that are new to the atmosphere
% ng_ch4_fuel_content  ... natural gas CH4 emissions (short tons/mmbtu)
% coal_ch4_fuel_content... coal CH4 emissions (short tons/mmbtu)
% eor_leakage          ... Leakage rate for CO2 stored in enhanced oil recovery

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Gather Data from ESC file
info = esc.gen_ann_adj;
info = filterInfo(info, esc);
info.value = info.value(:, ['Y', num2str(esc.year)]);



%natural gas CH4 emissions (short tons/mmbtu)
idx = strcmp(info.col_name{:,  'col_name'}, 'CH4_FUEL_CONTENT') ...
    & strcmp(info.genfuel{:, 'genfuel'}, 'ng');
assert(sum(idx) == 1, 'Error reading CH4 fuel content for natural gas plants')
eopt.ng_ch4_fuel_content = info.value{idx,:};

%coal CH4 emissions (short tons/mmbtu)
idx = strcmp(info.col_name{:, 'col_name'}, 'CH4_FUEL_CONTENT') ...
    & strcmp(info.genfuel{:, 'genfuel'}, 'coal');
assert(sum(idx) == 1, 'Error reading CH4 fuel content for coal plants')
eopt.coal_ch4_fuel_content = info.value{idx,:};

%global warming potential of methane (stons CO2/stons CH4)
idx = strcmp(info.col_name{:, 'col_name'}, 'CH4_WARMING_POTENTIAL');
assert(sum(idx) == 1, 'Error reading CH4 global warming potential')
eopt.ch4_gwp = info.value{idx,:};

% Leakage rate for CO2 stored in enhanced oil recovery
if strcmp(eopt.ccs, 'T') || strcmp(eopt.dac, 'T')
    idx = strcmp(info.col_name{:, 'col_name'}, 'EOR_LEAKAGE_RATE');
    assert(sum(idx) == 1, 'Error reading EOR CO2 leakage rate')
    eopt.eor_leakage = info.value{idx,:};
else
    eopt.eor_leakage = NaN;
end

% Fraction of biomass emissions CO2 emissions that count toward CO2e emissions.
idx = strcmp(info.col_name{:, 'col_name'}, 'BIOMASS_pctCO2e');
assert(sum(idx) == 1, 'Error reading biomass emissions multiplier')
eopt.bio_pctCO2e = info.value{idx,:};

%% Override ESC Inputs with Optional Input Data, if Present
p = inputParser; %create new structure p to store values

%set default values to those in the esc file
addOptional(p, 'ch4_gwp', eopt.ch4_gwp, @(x) isnumeric(x))
addOptional(p, 'bio_pctCO2e', eopt.bio_pctCO2e, @(x) isnumeric(x))
addOptional(p, 'ng_ch4_fuel_content', eopt.ng_ch4_fuel_content, @(x) isnumeric(x))
addOptional(p, 'coal_ch4_fuel_content', eopt.coal_ch4_fuel_content, @(x) isnumeric(x))
addOptional(p, 'eor_leakage', eopt.eor_leakage, @(x) isnumeric(x))

%override default values with user-specified inputs
parse(p, varargin{:})

%% Recalculate CO2e emission rates for various generators
% initialize
mpc.gen_aux{:, 'EMIS_CO2e'} = mpc.gen_aux{:, 'EMIS_CO2'};

%add methane emissions to ng and coal emissions
mpc.gen_aux{strcmp(mpc.genfuel, 'ng'), 'EMIS_CO2e'} = mpc.gen_aux{strcmp(mpc.genfuel, 'ng'), 'EMIS_CO2'} + (p.Results.ng_ch4_fuel_content * p.Results.ch4_gwp * mpc.gen_aux{strcmp(mpc.genfuel, 'ng'), 'HR'}); % (stons ch4/mmbtu) * (ston co2/ston ch4) * (mmbtu/mwh) => ston co2/mwh
mpc.gen_aux{strcmp(mpc.genfuel, 'coal'), 'EMIS_CO2e'} = mpc.gen_aux{strcmp(mpc.genfuel, 'coal'), 'EMIS_CO2'} + (p.Results.coal_ch4_fuel_content * p.Results.ch4_gwp * mpc.gen_aux{strcmp(mpc.genfuel, 'coal'), 'HR'}); % (stons ch4/mmbtu) * (ston co2/ston ch4) * (mmbtu/mwh) => ston co2/mwh
mpc.gen_aux{strcmp(mpc.genfuel, 'dac'), 'EMIS_CO2e'} = mpc.gen_aux{strcmp(mpc.genfuel, 'dac'), 'EMIS_CO2'} + (p.Results.ng_ch4_fuel_content * p.Results.ch4_gwp * mpc.gen_aux{strcmp(mpc.genfuel, 'dac'), 'HR'}); % (stons ch4/mmbtu) * (ston co2/ston ch4) * (mmbtu/mwh) => ston co2/mwh

%adjust biomass emissions for the fact that only a fraction of the CO2
%emitted is new to the atmosphere.
mpc.gen_aux{strcmp(mpc.genfuel, 'biomass'), 'EMIS_CO2e'} = p.Results.bio_pctCO2e * mpc.gen_aux{strcmp(mpc.genfuel, 'biomass'), 'EMIS_CO2'};

% For combined-heat-and-power (CHP) plants, assign only a fraction of
% carbon emissions to CO2e. This is becuase combined-heat-and-power plants
% are not pure power plants, but rather also produce heat for residential
% and industrial applications. Thus, not all CO2e emissions are for power generation.
mpc.gen_aux{:, 'EMIS_CO2e'} = mpc.gen_aux{:, 'EMIS_CO2e'} .* mpc.gen_aux{:, 'CHP_CO2_MULTI'};


% For CCS and DAC that store carbon using enhanced oil recovery (EOR),
% add the CO2 that is leaked back to the atmosphere.
if strcmp(eopt.ccs, 'T') || strcmp(eopt.dac, 'T')
    idx_eor = strcmp(mpc.gen_map{:, 'co2_storage_type'}, 'eor');
    mpc.gen_aux{idx_eor, 'EMIS_CO2e'} = mpc.gen_aux{idx_eor, 'EMIS_CO2e'} + p.Results.eor_leakage*mpc.gen_aux{idx_eor, 'STORAGE_CO2'}; %This works for dac as well, since CO2 coeff is positive but storage coeff is negative.
end
end

