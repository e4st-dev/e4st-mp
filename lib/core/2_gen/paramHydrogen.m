function [mpc, offer] = paramHydrogen(mpc, esc, offer, eopt)
%%PARAMHYDROGEN
% This file sets up hydrogen retrofits in E4ST. This file will update the
% mpc and offer matrices to allow every pre-existing natural gas combined
% cycle unit (ngcc) to be retrofited with a hydrogen turbine.
% 
% Retrofits are implemented in E4ST analagous to other buildable
% generators. The key difference is that they are linked to an existing
% generation unit through a capacity costraint. The capacity constraint
% specifies that the sum of the existing capacity and retrofited capacity
% must equal the existing capacity.
%
% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Paul Picciano (Resources for the Future), and Christoph Funke
% (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Column Numbers for MATPOWER Inputs
define_constants;

%% Define Cost Parameters for Units to Retrofit
if strcmp(eopt.hydrogen_retro, 'T')
    
    %identify units to retrofit
    idx_online = mpc.gen(:, GEN_STATUS) == 1;
    idx_retro = idx_online & ...
        (strcmp(mpc.gentype, 'ngcc') | strcmp(mpc.gentype, 'ngcc_new'))...
        & (mpc.newgen ~= 1); %include ngcc_new to allow newly built generators to be retrofited whent there is iteration 
    % & strcmp(mpc.gen_map{:, 'nation'}, 'us');
    
    %% Get data on existing units to be retrofitted
    [retrogen, retrogen.offer] = filter_e4st_gen(mpc, offer, idx_retro);
        
    %% Apply Retrofit Parameters
    %flag as newly built generators
    retrogen.newgen(:, :) = 1;
    
    %set different emission rates based on different levels of innovation
    emis_nox = 0.58;

    %set cost values for hydrogen combined cycle retrofits
    idx_retro_hcc = (strcmp(retrogen.gentype, 'ngcc') | strcmp(retrogen.gentype, 'ngcc_new'));
    retrogen.gentype(idx_retro_hcc) = {'hcc_ret'};
    retrogen.genfuel(idx_retro_hcc) = {'hydrogen'};
    retrogen.gen_aux{idx_retro_hcc, 'CAP_COST'} = 8.63/2; %half of natural gas turbine
    retrogen.gen_aux{idx_retro_hcc, 'FOM'} = 1.36; %same as ngccs
    retrogen.gen_aux{idx_retro_hcc, 'VOM'} = 2; %same as ngccs
    retrogen.gen_aux{idx_retro_hcc, 'HR'} = 6.4; % same as ngccs
    retrogen.gen_aux{idx_retro_hcc, 'EMIS_CO2'} = 0; 
    retrogen.gen_aux{idx_retro_hcc, 'EMIS_NOX'} = emis_nox; % 20 percent higher than ngccs
    retrogen.gen_aux{idx_retro_hcc, 'EMIS_SO2'} = 0;
%     retrogen.gen(idx_retro_hcc, PMAX) = retrogen.gen(idx_retro_hcc, PMAX); %leave unchanged
%     retrogen.offer(idx_retro_hcc, 2) = retrogen.offer(idx_retro_hcc, 2); %leave unchanged
    
    %% Add Retrofits to Main Data
    [mpc, offer] = append_e4st_gen(mpc, offer, retrogen);
    
    %% Update Offer and Gencost
    idx_retro = strcmp(mpc.gentype, 'hcc_ret');
    offer = updateOfferPrc(mpc, offer, idx_retro, eopt);
    mpc = updateGenCost(mpc, idx_retro, eopt);  

end

%% Apply Capacity Constraints:
% Retrofits
if strcmp(eopt.hydrogen_retro, 'T')
    idx_retro = ismember(mpc.gentype, {'hcc_ret'});
    ids = unique(mpc.gen_map{idx_retro, 'gen_id'});
    fprintf('Applying Hydrogen Retrofit Capacity Constraints (%d for %d total generators)\n', length(ids), sum(idx_retro))
    for id = ids'
        idx_id = mpc.gen_map{:, 'gen_id'} == id;
        
        if any(idx_id & idx_retro)
            iRetro = find(idx_id & idx_retro);
            iUnRetro = find(idx_id & ~idx_retro);
            cap = mpc.gen(iUnRetro(1), PMAX);
            map = double(idx_id'); %coefficients of generators in capacity constranit
            capMax = cap; % max capcity value
            capMin = -Inf; % min capcity value
            name = {['hydrogen_retrofit_', num2str(id)]};
            mpc = addCapLim(mpc, map, capMin, capMax, name);
        end
    end
end



%% Update MPC Dimensions
mpc = updateMPCDim(mpc);

end
