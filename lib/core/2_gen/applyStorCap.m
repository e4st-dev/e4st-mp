function [mpc] = applyStorCap(mpc)
% applyStorCap Applies capacity constraints to the CO2 storage steps 
% to which CCS and DAC plants can send CO2

% E4ST
% Copyright (c) 2012-2018 by Power System Engineering Research Center (PSERC)
% by Steven Witkin (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Apply Storage Quantity Constraints

stor_state_step = unique(mpc.gen_map(:, {'co2_storage_state', 'co2_storage_step'}));
stor_state_step(strcmp(stor_state_step{:, 'co2_storage_state'}, 'na'), :) = [];

stor_states = mpc.gen_map{:, 'co2_storage_state'};
stor_steps = mpc.gen_map{:, 'co2_storage_step'};
stor_step_tons = mpc.gen_map{:, 'co2_storage_step_cap'};
fprintf('Applying CO2 Storage Step Constraints (%d total)\n', height(stor_state_step))
for i = 1:height(stor_state_step)
    map = (strcmp(stor_state_step{i, 'co2_storage_state'}, stor_states) & strcmp(stor_state_step{i, 'co2_storage_step'}, stor_steps))';
    cap = unique(stor_step_tons(map')) / 8760; % Convert the cap from annual to hourly
    tocMin = -Inf;
    tocMax = cap;
    coeff = mpc.gen_aux{:, 'STORAGE_CO2'};
    name = strjoin([{'CO2_Storage_'}, stor_state_step{i, 'co2_storage_state'}, {'_'}, stor_state_step{i, 'co2_storage_step'}], '');
    mpc = addTOC(mpc, map, tocMin, tocMax, coeff, 1, name);
end

%% Update MPC Dimensions
mpc = updateMPCDim(mpc);

end