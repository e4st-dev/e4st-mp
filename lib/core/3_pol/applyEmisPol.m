function mpc = applyEmisPol(mpc, pol_val, pol, pol_type, idx_gen, pol_name, esc, eopt)
%applyEmisPol: Set emission caps and emissions prices for CO2, NOX, SO2, or PM25
% Also used to set caps and prices on CO2 storage

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Emission rate column in mpc.gen_aux
switch pol
    case 'CO2'
        emis_rate = 'EMIS_CO2';
    case 'CO2e'
        emis_rate = 'EMIS_CO2e';
    case 'CO2e_20y'
        emis_rate = 'EMIS_CO2e';
    case 'NOX'
        emis_rate = 'EMIS_NOX';
    case 'SO2'
        emis_rate = 'EMIS_SO2';
    case 'PM25'
        emis_rate = 'EMIS_PM25';
    case 'CO2_STOR'
        emis_rate = 'STORAGE_CO2';
    otherwise
        vfprintf(eopt.verbose, 'Error: Pollutant type not accepted\n')
end

%% Do not apply to DL
idx_gen = idx_gen & ~strcmp(mpc.genfuel, 'dl');

%% Apply Policy
if strcmp(pol, 'CO2e_20y')
    % Temp override of CO2e calculation using new value
    mpc2 = updateCO2e(mpc, esc, eopt, 'ch4_gwp', eopt.CH4_GWP_20y); %calculate emission rate with 20 year GWP of methane
    coeff = mpc2.gen_aux{:, emis_rate};
else
    coeff = mpc.gen_aux{:, emis_rate};
end

if strcmp(pol_type, 'cap')
    cap = pol_val / 8760; % Convert the cap from annual to hourly
    mpc = addTOC(mpc, idx_gen', -Inf, cap, coeff, 1, {pol_name});
    if strcmp(pol, 'CO2_STOR')
        vfprintf(eopt.verbose, 'Annual CO2 storage cap of %.2f applied to %d generators\n', pol_val, sum(idx_gen));
    else
        vfprintf(eopt.verbose, 'Annual %s emissions cap of %.2f applied to %d generators\n', pol, pol_val, sum(idx_gen));
    end

elseif strcmp(pol_type, 'prc') || strcmp(pol_type, 'price')
    emis_prc = pol_val;
    mpc.gen_aux{idx_gen, 'PTC'} = mpc.gen_aux{idx_gen, 'PTC'} - (coeff(idx_gen) * emis_prc);
    mpc = updateGenCost(mpc, idx_gen, eopt);
    if strcmp(pol, 'CO2_STOR')
        vfprintf(eopt.verbose, 'CO2 storage tax of $%f applied to %d generators\n', pol_val, sum(idx_gen));
    else
        vfprintf(eopt.verbose, '%s emission price (i.e. tax) of $%f applied to %d generators\n', pol, pol_val, sum(idx_gen));
    end
end
