function [mpc, offer] = retireGen(mpc, offer, esc, result, eopt)
%retireGen: Retire generators that are smaller than theta or not used

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Initialize
define_constants;
idx_online = mpc.gen(:, GEN_STATUS) == 1;

%% Select generators to retire
idx_can_ret = idx_online & ~strcmp(mpc.genfuel, 'dl'); % Only retire online generators, and not dispatchable load
idx_can_ret_nodac = idx_online & ~ismember(mpc.genfuel, {'dl', 'dac'});
idx_can_ret_dac = idx_online & strcmp(mpc.genfuel, 'dac');

%% Endogenous Retirements
if esc.year >= esc.first_ret_yr && ~isempty(result)
    cap_used = result.reserve.qty.Rp_pos; % Capacity of online generators
    % In e4st_solve, short term storage units are made to have equal Rp_neg
    % and Rp_pos, so this works fine.
    idx_negpmin = mpc.gen(:, PMIN) < 0; % & mpc.gen(:, PMAX) <= 0;
    cap_used(idx_negpmin) = result.reserve.qty.Rp_neg(idx_negpmin);
    cap_shutdown = sum(offer(idx_online, 2)-cap_used(idx_online, 1)); % check if offline generators included here, also not exact for DAC

    % Baseload Coal/Hydro Retirements
    %if any(strcmp(mpc.gen_map.Properties.VariableNames, 'op_type')) && any(strcmp(mpc.gen_map{:, 'op_type'}, 'baseload')) && strcmp(esc.ret_bl, 'T')
    %    idx_bl = strcmp(mpc.gen_map{:, 'op_type'}, 'baseload');
    %    idx_dp = strcmp(mpc.gen_map{:, 'op_type'}, 'dispatchable');
    %    idx_check = offer(:, 2) > 0;
    %    ids = unique(mpc.gen_map{idx_dp & idx_check, 'gen_id'});
    %    for id = ids'
    %        idx_id = mpc.gen_map{:, 'gen_id'} == id;
    %        if sum(idx_id) ~= 2
    %            continue
    %        end
    %        ratio = cap_used(idx_id & idx_dp, 1) / offer(idx_id & idx_dp, 2);
    %        cap_used(idx_id & idx_bl) = cap_used(idx_id & idx_bl) * ratio;
    %    end
    %end

    % Retirement: Set PositiveReserveCap to used capacity
    offer(idx_can_ret_nodac, 2) = cap_used(idx_can_ret_nodac);
    mpc.gen(idx_can_ret_nodac, PMAX) = cap_used(idx_can_ret_nodac);
    mpc.gen(idx_can_ret_nodac, PG) = cap_used(idx_can_ret_nodac);
    
    % For DAC
    offer(idx_can_ret_dac, 4) = cap_used(idx_can_ret_dac);
    mpc.gen(idx_can_ret_dac, PMIN) = -cap_used(idx_can_ret_dac);
    mpc.gen(idx_can_ret_dac, PG) = -cap_used(idx_can_ret_dac);
    

    % Switch off generators with zero capacity (including those with capacity less than theta)
    % Storage units use the retire_theta, DAC uses a smaller retirement criteria
    idx_ret_nodac = idx_can_ret_nodac & (cap_used < esc.retire_theta);
    idx_ret_dac = idx_can_ret_dac & (cap_used < esc.retire_theta/10);
    idx_retire = idx_ret_nodac | idx_ret_dac;
    mpc.gen(idx_retire, GEN_STATUS) = 0;

    if isfield(mpc, 'gen_desc')
        mpc.gen_desc{idx_retire, 'RET_REASON'} = {'endogenous'};
        mpc.gen_desc{idx_retire, 'RET_YR'} = esc.year;
        mpc.gen_desc{idx_retire, 'RET_SIMYR'} = {eopt.sim_yr};
    end

    vfprintf(eopt.verbose, '%d generators, totalling %.1f MW, have been endogenously shutdown and retired\n', sum(idx_retire), cap_shutdown);
end

%% Exogenous Retirements
if isfield(mpc, 'gen_desc') && any(idx_can_ret & mpc.gen_desc{:, 'P_RET_YR'} <= esc.year)
    idx_retire = idx_can_ret & mpc.gen_desc{:, 'P_RET_YR'} <= esc.year;
    cap_shutdown = sum(mpc.gen(idx_retire, PMAX));

    mpc.gen(idx_retire, GEN_STATUS) = 0;
    %offer(idx_retire, 2) = 0;
    %mpc.gen(idx_retire, PMAX) = 0;
    %mpc.gen(idx_retire, PG) = 0;

    if isfield(mpc, 'gen_desc')
        mpc.gen_desc{idx_retire, 'RET_REASON'} = {'exogenous'};
        %mpc.gen_desc{idx_retire, 'RET_YR'} = esc.year;
        mpc.gen_desc{idx_retire, 'RET_YR'} = mpc.gen_desc{idx_retire, 'P_RET_YR'};
        mpc.gen_desc{idx_retire, 'RET_SIMYR'} = {eopt.sim_yr};
    end

    vfprintf(eopt.verbose, '%d generators, totalling %.1f MW, have been exogenously shutdown and retired\n', sum(idx_retire), cap_shutdown);
end
