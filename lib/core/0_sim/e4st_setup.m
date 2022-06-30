function [mpc, esc, offer, contab, eopt] = e4st_setup(mpc, esc, eopt)
% e4st_setup prepares e4st_simulation

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% ----------DISPLAY KEY SIMULATION INFO---------
if eopt.verbose.verbose
    fprintf('\n---------------------------------------------------------\n')
    fprintf('Preparing E4ST Simulation:\n')
    if eopt.iter.status == 1
        fprintf('\tSim Name: %s\n', eopt.sim_name_yr_iter)
    else
        fprintf('\tSim Name: %s\n', eopt.sim_name_yr)
    end
    fprintf('\tCase Name: %s\n', eopt.case_name)
    fprintf('\tGrid Name: %s\n', eopt.grid_name)
    fprintf('\tYear: %d\n', esc.year);
    if strcmp(eopt.prl, 'T')
        fprintf('\tLoad Treatment: Price Responsive Load\n')
    elseif strcmp(eopt.default, 'T')
        fprintf('\tLoad Treatment: Fixed Load (Default Simulation)\n')
    else
        fprintf('\tLoad Treatment: Fixed Load (Result Simulation)\n')
    end
    if isfield(eopt, 'mpc_file')
        fprintf('\tMPC file: %s\n', eopt.mpc_file)
    end
    if isfield(eopt, 'esc_file')
        fprintf('\tESC file: %s\n', eopt.esc_file)
    end
    fprintf('---------------------------------------------------------\n\n')
end

%% -------------- 1. SETUP: GRID --------------
[mpc, esc] = setup_grid(mpc, esc, eopt);

%% -------------- 2. SETUP: GENERATORS --------------
[mpc, offer, eopt] = setup_gen(mpc, esc, eopt);

%% -------------- 3. SETUP: POLICY --------------
[mpc, offer] = setup_policy(mpc, offer, esc, eopt);

%% -------------- 4. SETUP: LOAD --------------
[mpc, contab] = setup_load(mpc, esc, eopt, 1);

%% Sync Offer & Gencost
offer = syncOffer(mpc, offer, eopt);

%% Override Offer & Gencost for Simulation
if isfield(esc, 'gen_ann_adj')
    [mpc, offer] = applyAnnAdj(mpc, offer, esc, 'offer', eopt); % Quasi-Fixed Costs
    [mpc, offer] = applyAnnAdj(mpc, offer, esc, 'gencost', eopt); % Short-Run Marginal Costs
end

