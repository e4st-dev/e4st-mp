function [mpc, contab] = setup_load(mpc, esc, eopt, make_cont)
% setup_load prepare electricity demand values for simulation

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

% When setup_load is run in setup_policy to get load amounts, contab is not
% needed. Save time by skipping it using make_cont
vfprintf(eopt.verbose, '\n\n-----Running E4ST Set-Up Module: Demand (Load)-----\n\n')

%% Scale base loads by annual growth rate
if isfield(esc, 'load_growth') && strcmp(eopt.update_gens, 'T')
    mpc = scaleBaseLoad(mpc, esc, eopt);
end

%% Scale load in contingency hours
if isfield(esc, 'load_shape') && strcmp(eopt.update_gens, 'T')
    mpc = shapeLoad(mpc, esc, eopt);
end

%% Scale base loads to match true values
if isfield(esc, 'load_match') && strcmp(eopt.update_gens, 'T')
    mpc = matchLoad(mpc, esc, eopt);
end

%% Add to hourly loads, as necessary
if isfield(esc, 'load_add') && strcmp(eopt.update_gens, 'T')
    [mpc] = addLoad(mpc, esc, eopt);     
end

%%
if make_cont
    %% create MATPOWER contingency table based on hourly load
    [mpc, contab] = create_contab(mpc, esc, eopt);
    
    
    %% Price responsive load
    if strcmp(eopt.prl, 'T') && ~strcmp(eopt.default, 'T')
        [mpc, contab] = makePrlNew(mpc, contab, esc, eopt);
    end
else
    contab = [];
end
