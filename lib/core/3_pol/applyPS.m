function mpc = applyPS(mpc, esc, ps_val, ps_type, cred_type, mapGen, mapLoad, name, eopt)
%applyPS: Set portfolio standard constraints (RPS/CES)

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).


%NOTES: mapGen = idx_gen and map_load = idx_dl

%% DAC does not earn RPS credits
if strcmp(cred_type, 'rps')
    mapGen = mapGen & ~strcmp(mpc.genfuel, 'dac'); % if dac doesn't explicitly get credits, assume it doesn't earn and just has to buy
end

%% Calculate Credit
if strcmp(cred_type, 'ces')
    gen_er = mpc.gen_aux{:, 'EMIS_CO2e'};
    credit = 1 - (gen_er / bm_er);
    credit(strcmp(mpc.genfuel, 'biomass')) = min(.5*(2204 / 2000), credit(strcmp(mpc.genfuel, 'biomass')));
    %credit(strcmp(mpc.genfuel,'biomass')) = .5 * (2204 / 2000);
    credit(credit < 0) = 0;
    
elseif strcmp(cred_type, 'ces_coal')
    bm_er = .82 * 2204 / 2000; % Coal Emis Rate (Short Tons CO2e)
    gen_er = mpc.gen_aux{:, 'EMIS_CO2e'};
    credit = 1 - (gen_er / bm_er);
    credit(strcmp(mpc.genfuel, 'biomass')) = min(.5*(2204 / 2000), credit(strcmp(mpc.genfuel, 'biomass')));
    credit(strcmp(mpc.gentype, 'coal_cofire')) = .5 * .15;
    credit(credit < 0) = 0;
    credit(strcmp(mpc.genfuel, 'dac')) = - gen_er(strcmp(mpc.genfuel, 'dac')) / bm_er;

elseif strcmp(cred_type, 'clean_futures_act_2021')
    if esc.year <= 2030
        bm_er = 0.82 * 2204 / 2000;
    elseif esc.year <= 2031
        bm_er = 0.736 * 2204 / 2000;
    elseif esc.year <= 2032
        bm_er = 0.652 * 2204 / 2000;
    elseif esc.year <= 2033
        bm_er = 0.568 * 2204 / 2000;
    elseif esc.year <= 2034
        bm_er = 0.484 * 2204 / 2000;
    else
        bm_er = 0.4 * 2204 / 2000;        
    end
    mpc2 = updateCO2e(mpc, esc, eopt, 'ch4_gwp', 86);%calculate emission rate with 20 year GWP of methane
    gen_er = mpc2.gen_aux{:, 'EMIS_CO2e'};
    credit = 1 - (gen_er / bm_er);
    credit(credit < 0) = 0;
    credit(strcmp(mpc.genfuel, 'dac')) = - gen_er(strcmp(mpc.genfuel, 'dac')) / bm_er;
    
    
elseif strcmp(cred_type, 'ces_ngcc')
    bm_er = .4 * 2204 / 2000; % New NGCC Emis Rate (Short Tons CO2e)
    gen_er = mpc.gen_aux{:, 'EMIS_CO2e'};
    credit = 1 - (gen_er / bm_er);
    credit(strcmp(mpc.genfuel, 'biomass')) = min(.5*(2204 / 2000), credit(strcmp(mpc.genfuel, 'biomass')));
    credit(strcmp(mpc.gentype, 'coal_cofire')) = .5 * .15;
    credit(credit < 0) = 0;
    credit(strcmp(mpc.genfuel, 'dac')) = - gen_er(strcmp(mpc.genfuel, 'dac')) / bm_er;
    
elseif strcmp(cred_type, 'rps')
    % default is NG-CCS gets <1 credit according to CO2e
    % Does nothing if ng-ccs not included in policy in esc file.
    credit = ones(size(mapGen, 1), 1);
    idx_ngccs = strcmp(mpc.gentype, 'ngccccs_new');
    credit(idx_ngccs) = 1 - mpc.gen_aux{idx_ngccs, 'EMIS_CO2e'} / (.6 * 2204 / 2000);
    credit(credit < 0) = 0;
else
    credit = ones(size(mapGen, 1), 1);
end

%% Apply portfolio standard
if strcmp(ps_type, 'perc')
    mapPS = mapGen | (mapLoad ~= 0);
    credit = credit .* mapGen;
    coeff = (ps_val * -mapLoad) - credit; % load * ps_req <= gen (total_output <= K). There is a negative sign before dl.
    cap = 0;
    mpc = addTOC(mpc, mapPS', -Inf, cap, coeff, 1, {name});
    vfprintf(eopt.verbose, 'PS percentage requirement of %.3f applied to %d generators and %d electricity demand units\n', mean(mapLoad(mapLoad > 0)), sum(mapGen), sum(mapLoad ~= 0));

elseif strcmp(ps_type, 'qty')
    assert(~strcmp(eopt.dac, 'T'), "Rps/ces qty policies not set up for use with DAC yet");
    mapPS = mapGen;
    coeff = credit .* mapGen;
    cap = ps_val / 8760; % hourly output requirement

    define_constants;
    max_output = sum(mpc.gen(mapGen, PMAX).*(mpc.availability_factor(mapGen, :) * esc.hrs_map{:, 'hours'})) / 8760;
    if max_output >= cap % Check if feasible
        mpc = addTOC(mpc, mapPS', cap, Inf, coeff, 1, {name});
        vfprintf(eopt.verbose, 'PS quantity requirement of %.3f applied to %d generators\n', ps_val, sum(mapGen));
    else
        vfprintf(eopt.verbose, 'Warning: Not enough generator capacity to meet portfolio standard\n')
    end

elseif strcmp(ps_type, 'qty_max')
    assert(~strcmp(eopt.dac, 'T'), "Rps/ces qty_max policies not set up for use with DAC yet");
    mapPS = mapGen;
    coeff = credit .* mapGen;
    cap = ps_val / 8760; % hourly output requirement

    define_constants;
    mpc = addTOC(mpc, mapPS', -Inf, cap, coeff, 1, {name});
    vfprintf(eopt.verbose, 'PS quantity requirement of %.3f applied to %d generators\n', ps_val, sum(mapGen));

elseif strcmp(ps_type, 'prc') || strcmp(ps_type, 'price')
    % Storage and DAC have to pay at 100% for credits
    credit(strcmp(mpc.genfuel, 'dac')) = credit(strcmp(mpc.genfuel, 'dac')) + 1;
    credit(strcmp(mpc.genfuel, 'storage')) = -0.15/0.85;
    
    ps_sub = ones(length(mapGen), 1) .* mapGen .* credit * ps_val;
    mpc.gen_aux{mapGen, 'PTC'} = mpc.gen_aux{mapGen, 'PTC'} + ps_sub(mapGen);
    mpc = updateGenCost(mpc, mapGen, eopt);
    vfprintf(eopt.verbose, 'PS price of %.3f applied to %d generators\n', ps_val, sum(mapGen));
end
