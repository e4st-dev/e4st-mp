function [mpc, offer] = matchCap(mpc, esc, offer, eopt)
%matchCap: Set capacity constraints

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Prepare Data
info = esc.cap_match;
info = filterInfo(info, esc);

info.value = info.value(:, ['Y', num2str(esc.year)]);

if isfield(info, 'joint')
    info_list = info.joint{:, :};
else
    info_list = (1:height(info.value))';
end

%% Match Capacity
define_constants;
idx_online = mpc.gen(:, GEN_STATUS) == 1;
cap_num = 0;

for i = unique(info_list)'
    idx_info = i == info_list;
    cur_info = filterStruct(info, idx_info);

    idx_gen = getInfoIdx(mpc, esc, cur_info, 'gen');
    idx_gen = idx_gen & idx_online;

    if ~any(idx_gen)
        continue
    end

    method = unique(cur_info.method{:, :});
    value = unique(cur_info.value{:, :});
    cap_num = cap_num + 1;

    switch char(method)
        case 'scale'
            if sum(mpc.gen(idx_gen, PMAX)) > 0
                scaler = value / sum(mpc.gen(idx_gen, PMAX));
                mpc.gen(idx_gen, PMAX) = mpc.gen(idx_gen, PMAX) * scaler;
                mpc.gen(idx_gen, PG) = mpc.gen(idx_gen, PG) * scaler;
                mpc.gen(idx_gen, PMIN) = mpc.gen(idx_gen, PMIN) * scaler;
                offer(idx_gen, 2) = mpc.gen(idx_gen, PMAX);
                
                % if capacity scaler is zero, count as exogenous retirement
                if scaler == 0 && isfield(mpc, 'gen_desc')
                    mpc.gen_desc{idx_gen, 'RET_REASON'} = {'exogenous'};
                    mpc.gen_desc{idx_gen, 'RET_YR'} = esc.year;
                    mpc.gen_desc{idx_gen, 'RET_SIMYR'} = {eopt.sim_yr};
                    mpc.gen(idx_gen, GEN_STATUS) = 0;
                end
            else
                vfprintf(eopt.verbose, 'Capacity not matched by scaling because no existing capacity\n')
            end

        case 'retire'
            mpc.gen(idx_gen, GEN_STATUS) = 0;

            if isfield(mpc, 'gen_desc')
                mpc.gen_desc{idx_gen, 'RET_REASON'} = {'exogenous'};
                mpc.gen_desc{idx_gen, 'RET_YR'} = esc.year;
                mpc.gen_desc{idx_gen, 'RET_SIMYR'} = {eopt.sim_yr};
            end

        case 'max_perc'
            cap_cur = sum(mpc.gen(idx_gen, PMAX));
            cap_match = cap_cur * value;
            idx_gen = idx_gen & ~strcmp(mpc.genfuel, 'dl');
            min = -Inf;
            max = cap_match;
            cap_name = {['CAPMATCH_', num2str(i), '_', char(unique(cur_info.area{:, 'area'})), '_', char(unique(cur_info.area{:, 'subarea'}))]};
            mpc = addCapLim(mpc, idx_gen', min, max, cap_name);

        case 'max_qty'
            cap_match = value;
            idx_gen = idx_gen & ~strcmp(mpc.genfuel, 'dl');
            min = -Inf;
            max = cap_match;
            cap_name = {['CAPMATCH_', num2str(i), '_', char(unique(cur_info.area{:, 'area'})), '_', char(unique(cur_info.area{:, 'subarea'}))]};
            mpc = addCapLim(mpc, idx_gen', min, max, cap_name);

        case 'min_perc'
            cap_cur = sum(mpc.gen(idx_gen, PMAX));
            cap_match = cap_cur * value;
            idx_gen = idx_gen & ~strcmp(mpc.genfuel, 'dl');
            min = cap_match;
            max = Inf;
            cap_name = {['CAPMATCH_', num2str(i), '_', char(unique(cur_info.area{:, 'area'})), '_', char(unique(cur_info.area{:, 'subarea'}))]};
            mpc = addCapLim(mpc, idx_gen', min, max, cap_name);

        case 'min_qty'
            cap_match = value;
            idx_gen = idx_gen & ~strcmp(mpc.genfuel, 'dl');
            min = cap_match;
            max = Inf;
            cap_name = {['CAPMATCH_', num2str(i), '_', char(unique(cur_info.area{:, 'area'})), '_', char(unique(cur_info.area{:, 'subarea'}))]};
            mpc = addCapLim(mpc, idx_gen', min, max, cap_name);
    end
end
end
