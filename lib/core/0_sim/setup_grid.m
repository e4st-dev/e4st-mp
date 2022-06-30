function [mpc, esc] = setup_grid(mpc, esc, eopt)
% setup_grid prepare grid data for simulation

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).


vfprintf(eopt.verbose, '\n\n-----Running E4ST Set-Up Module: Grid-----\n\n')

%% Extract and remove isolated buses
if isfield(esc, 'island')
    [mpc, esc] = selIsland(mpc, esc, esc.island, eopt); % 0 = Remove all islands
end

%% Increase all Line limits and Interface limits
if isfield(eopt, 'trans_increase')
    define_constants;
    factor = eopt.trans_increase{:, ['Y' num2str(esc.year)]};
    mpc.branch(:, RATE_A) = mpc.branch(:, RATE_A) * factor;
    mpc.branch(:, RATE_B) = mpc.branch(:, RATE_B) * factor;
    mpc.branch(:, RATE_C) = mpc.branch(:, RATE_C) * factor;
    vfprintf(eopt.verbose, "All branch line limits increased by a factor of %f \n", factor);
    
    if isfield(esc, 'interface_lim')
        esc.interface_lim.value{:, :} = esc.interface_lim.value{:, :} * factor;
        vfprintf(eopt.verbose, "All interface limits increased by a factor of %f \n", factor);
    end
else
    vfprintf(eopt.verbose, "Branch line and interface limits not increased");
end

%% Add interface constraints
if isfield(esc, 'interface_lim')
    mpc.if = [];
    mpc = addInterface(mpc, esc, eopt);
end

%% Active DC Lines
if isfield(mpc, 'dcline') && ~isempty(mpc.dcline)
    if toggle_dcline(mpc, 'status') %remove any existing DC-line structures
        mpc = toggle_dcline(mpc, 'off');
    end
    mpc = toggle_dcline(mpc, 'on'); %re-add the structures
end
