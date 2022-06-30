function [mpc, offer] = paramDAC(mpc, esc, offer, eopt)
% paramDAC set up direct air capture "generators"

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Steven Witkin (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%%
define_constants;

% Storage options per plant
stor_steps = esc.stor_steps;
n_stor = eopt.n_stor;

%% CO2 Storage Options

idx_dac = ismember(mpc.genfuel, {'dac'});
% We don't need to re-create storage options for units that already have
% them from previous iteration or anything else. Assuming default is STEP0
idx_dac = idx_dac & strcmp(mpc.gen_map{:, 'co2_storage_step'}, 'STEP0'); %only the unassigned DAC

gen_states = mpc.gen_map{:, 'state'};
prod_states = unique(gen_states(idx_dac));
for prod_state = prod_states'
    idx_gen = idx_dac & strcmp(gen_states, prod_state); %gens in the state
    gen_stor_steps = stor_steps(strcmp(stor_steps{:, 'ProdState'}, prod_state), :);
    gen_stor_steps = sortrows(gen_stor_steps, 'total_cost_2013'); %order from cheapest to most expensive
    gen_stor_steps = gen_stor_steps(1:n_stor, :); %only keep the first n 
    [mpc_subset, mpc_subset.offer] = filter_e4st_gen(mpc, offer, idx_gen);
    stor_states = gen_stor_steps{:, 'StorState'}; % storage state options for the gens in prod_state
    fprintf('Creating DAC CO2 Storage Options by Production State: %s (%d generators x %d storage states = %d total)\n', ...
        char(prod_state), sum(idx_gen), length(stor_states), sum(idx_gen)*length(stor_states))
    for i_stor_state = 1:length(stor_states)
        stor_state = stor_states(i_stor_state);
        %fprintf('\tStorage State: %s\n',char(stor_state))
        if i_stor_state == 1
            mpc.gen_map{idx_gen, 'co2_storage_state'} = stor_state;
            mpc.gen_map{idx_gen, 'co2_storage_step'} = gen_stor_steps{i_stor_state, 'StepName'};
            mpc.gen_map{idx_gen, 'co2_storage_step_cap'} = gen_stor_steps{i_stor_state, 'step_tons'};
            mpc.gen_map{idx_gen, 'co2_storage_cost'} = gen_stor_steps{i_stor_state, 'total_cost_2013'} * mpc.gen_aux{idx_gen, 'STORAGE_CO2'};
            mpc.gen_map{idx_gen, 'co2_storage_type'} = gen_stor_steps{i_stor_state, 'stor_type'};
        else
            mpc_subset.gen_map{:, 'co2_storage_state'} = stor_state;
            mpc_subset.gen_map{:, 'co2_storage_step'} = gen_stor_steps{i_stor_state, 'StepName'};
            mpc_subset.gen_map{:, 'co2_storage_step_cap'} = gen_stor_steps{i_stor_state, 'step_tons'};
            mpc_subset.gen_map{:, 'co2_storage_cost'} = gen_stor_steps{i_stor_state, 'total_cost_2013'} * mpc.gen_aux{idx_gen, 'STORAGE_CO2'};
            mpc_subset.gen_map{:, 'co2_storage_type'} = gen_stor_steps{i_stor_state, 'stor_type'};
            [mpc, offer] = append_e4st_gen(mpc, offer, mpc_subset);
        end
    end
end

% DAC units "generate" negative electricity, so their gencost is actually
% negative, since multiplying negative gen by negative gencost gives
% positive cost. So, their additional cost will actually be added to PTC so
% that additional cost can LOWER gencost, increasing it in magnitude.
% However, the negative storage value already does the negation, so we
% subtract like normal....
idx_dac = ismember(mpc.genfuel, {'dac'});
% could be replaced by applyPTC function call, but I think inline is quicker.
% PTC gets reset in setup_gen, so apply to all DAC
mpc.gen_aux{idx_dac, 'PTC'} = mpc.gen_aux{idx_dac, 'PTC'} - mpc.gen_map{idx_dac, 'co2_storage_cost'};
mpc = updateGenCost(mpc, idx_dac, eopt);

%% Apply Capacity Constraints:
% Not applying capacity constraints on units right now since all are new. 
% If need to, examples are in paramCCS.

%% Apply Storage Quantity Constraints
% implemented in applyStorCap after both CCS and DAC are setup

%% Update MPC Dimensions
mpc = updateMPCDim(mpc);

end