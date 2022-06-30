function mpc = addInterface(mpc, esc, eopt)
% addInterface Add interface limit (iflim) power flow constraint across set of specified branches

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Paul Picciano (Resources for the Future) and Biao Mao (Rensselaer Polytechnic Institute)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Initialize
info = esc.interface_lim;
info = filterInfo(info, esc);

if isfield(info, 'joint')
    info_list = info.joint{:, :};
else
    info_list = (1:height(info.value));
end

[if_lims, if_lims_hrs, if_map, if_name] = deal([]);
define_constants;

%% Prepare iflim parameters
if_num = 0;
for i = unique(info_list)'
    idx_info = i == info_list;

    cur_info = filterStruct(info, idx_info);

    % Idx for branches with 'f_bus' in 'from' area
    cur_info.area = cur_info.f_area;
    cur_info.area.Properties.VariableNames = {'area', 'subarea'};
    idx_br_f1 = getInfoIdx(mpc, esc, cur_info, 'branch');

    % Idx for branches with 't_bus' in 'from' area
    tmp_mpc = mpc;
    tmp_mpc.branch(:, 1) = tmp_mpc.branch(:, 2);
    idx_br_f2 = getInfoIdx(tmp_mpc, esc, cur_info, 'branch');

    % Idx for branches with 'f_bus' in 'to' area
    cur_info.area = cur_info.t_area;
    cur_info.area.Properties.VariableNames = {'area', 'subarea'};
    idx_br_t2 = getInfoIdx(mpc, esc, cur_info, 'branch');

    % Idx for branches with 't_bus' in 'to' area
    tmp_mpc = mpc;
    tmp_mpc.branch(:, 1) = tmp_mpc.branch(:, 2);
    idx_br_t1 = getInfoIdx(tmp_mpc, esc, cur_info, 'branch');

    % Branches with f_bus in "from" area and t_bus in "to" area
    idx_1 = idx_br_f1 & idx_br_t1 & mpc.branch(:, BR_STATUS);

    % Branches with t_bus in "from" area and f_bus in "to" area
    idx_2 = idx_br_f2 & idx_br_t2 & mpc.branch(:, BR_STATUS);

    % Check if any branches identified for interface flow limit
    if ~any(idx_1 | idx_2)
        continue
    end

    if_num = if_num + 1;

    % Specified flow limits are multiplied by negative one to reflect reverse direction
    idx_br = [find(idx_1); -1 * find(idx_2)];

    nh = size(esc.hrs_map, 1);
    direction = unique(cur_info.direction{:, :});

    for i = 1:nh
        lim = unique(cur_info.value{:, i});

        % Direction = 1 means limit applies one direction ("from to "to" area)
        % Direction = 2 means limit applies both directions ("from" to "to" area; "to" to "from" area)
        if lim < 0 && direction == 1
            % Min/max multipliers switched if specified limit is negative one direction (otherwise infeasible, I think)
            max_lim = lim * unique(cur_info.min{:, :});
            min_lim = lim * unique(cur_info.max{:, :});
        elseif lim < 0 && direction == 2
            % Min/max multipliers set to infinity if specified limit is negative (otherwise infeasible, I think)
            max_lim = Inf;
            min_lim = -Inf;
        else
            max_lim = lim * unique(cur_info.max{:, :});
            min_lim = lim * unique(cur_info.min{:, :});
        end

        % Set upper/lower bounds to infinity if user-specified limit is large
        if abs(max_lim) >= 9999
            max_lim = max_lim * Inf;
        end
        if abs(min_lim) >= 9999
            min_lim = min_lim * Inf;
        end

        if_lims_hrs(:, :, i) = [if_num, min_lim, max_lim];
    end

    % Store paramters for iflim constraint
    if_lims = [if_lims; if_lims_hrs];
    if_map = [if_map; [repmat(if_num, size(idx_br, 1), 1), idx_br]];
    if_name = [if_name; unique(cur_info.name{:, 'name'})];

end

%% Update iflim parameters
% Initialize
if ~isempty(if_map)
    if ~isfield(mpc, 'if')
        mpc.if = [];
    end
    if isempty(mpc.if)
        mpc.if.map = [];
        mpc.if.lims = [];
        mpc.if.name = [];
    end
end

% Apply
if ~isempty(if_map) && ~isempty(if_lims)
    if toggle_iflims(mpc, 'status') % Reset
        mpc = toggle_iflims(mpc, 'off');
    end
    mpc.if.map = [mpc.if.map; if_map];
    mpc.if.lims = [mpc.if.lims; if_lims];
    mpc.if.name = [mpc.if.name; if_name];
    mpc = toggle_iflims(mpc, 'on');
end

%% Display progress information
vfprintf(eopt.verbose, '%d interface constraint(s) added across %d branches\n', size(if_lims, 1), size(if_map, 1))
