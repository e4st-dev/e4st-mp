function res_constraint = result_constraints(res, mpc, esc)
% result_constraints Summarize results of E4ST constraints

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Results
define_constants;

%% Total Output Constraints
if isfield(res, 'total_output')
    toc_res = array2table(sum(sum(res.total_output.map > 0, 2), 3), 'VariableNames', {'map_count'});
    toc_res{:, 'min'} = res.total_output.min;
    toc_res{:, 'max'} = res.total_output.max;
    toc_res{:, 'type'} = res.total_output.type;
    toc_res{:, 'qty'} = res.total_output.qty;
    toc_res{:, 'mu'} = res.total_output.mu;
    toc_res{:, 'check_min'} = toc_res{:, 'mu'} .* (toc_res{:, 'qty'} - toc_res{:, 'min'});
    toc_res{:, 'check_max'} = toc_res{:, 'mu'} .* (toc_res{:, 'max'} - toc_res{:, 'qty'});

    if isfield(mpc.total_output, 'name')
        toc_res = [array2table(mpc.total_output.name, 'VariableNames', {'name'}), toc_res];
    end
    res_constraint.toc = toc_res;
else
    res_constraint.toc = [];
end

%% Capacity Limit Constraints
if isfield(res, 'caplim')
    caplim_res = array2table(sum(res.caplim.map > 0, 2), 'VariableNames', {'map_count'});
    caplim_res{:, 'min'} = res.caplim.min;
    caplim_res{:, 'max'} = res.caplim.max;
    caplim_res{:, 'qty'} = res.caplim.qty;
    caplim_res{:, 'mu'} = res.caplim.mu;
    caplim_res{:, 'check_min'} = caplim_res{:, 'mu'} .* (caplim_res{:, 'qty'} - caplim_res{:, 'min'});
    caplim_res{:, 'check_max'} = caplim_res{:, 'mu'} .* (caplim_res{:, 'max'} - caplim_res{:, 'qty'});

    if isfield(mpc.caplim, 'name')
        caplim_res = [array2table(mpc.caplim.name, 'VariableNames', {'name'}), caplim_res];
    end
    res_constraint.caplim = caplim_res;
end

%% Interface Limit Constraints
if isfield(res.base, 'if')
    nh = size(esc.hrs_map, 1);
    nc = max(res.base.if.map(:, 1));
    iflim_res = [];

    % hourly results for each interface
    for i = 1:nh
        if i == 1
            tmp_res = res.base;
        else
            tmp_res = res.cont(i-1);
        end

        tmp_iflim_res = array2table((1:nc)', 'VariableNames', {'if_num'});
        tmp_iflim_res{:, 'hr_num'} = i * ones(nc, 1);
        map_count = tabulate(tmp_res.if.map(:, 1));
        tmp_iflim_res{:, 'map_count'} = map_count(:, 2);
        tmp_iflim_res{:, 'min'} = tmp_res.if.lims(:, 2);
        tmp_iflim_res{:, 'max'} = tmp_res.if.lims(:, 3);
        tmp_iflim_res{:, 'qty'} = tmp_res.if.P;
        tmp_iflim_res{:, 'mu_min'} = tmp_res.if.mu.l / esc.hrs_map{i, 'probability'};
        tmp_iflim_res{:, 'mu_max'} = tmp_res.if.mu.u / esc.hrs_map{i, 'probability'};
        tmp_iflim_res{:, 'mu'} = max([tmp_iflim_res{:, 'mu_min'}, tmp_iflim_res{:, 'mu_max'}], [], 2);
        tmp_iflim_res{:, 'flow_mwh'} = tmp_iflim_res{:, 'qty'} * esc.hrs_map{i, 'hours'};
        for j = unique(tmp_res.if.map(:, 1))'
            idx = tmp_res.if.map(:, 1) == j;
            tmp_iflim_res{j, 'merch_surplus'} = sum(res.user_result.branch_res.ms_ac_hourly(idx, i)) * esc.hrs_map{i, 'hours'};
        end
        tmp_iflim_res{:, 'check_min'} = tmp_iflim_res{:, 'mu_min'} .* (tmp_iflim_res{:, 'min'} - tmp_iflim_res{:, 'qty'});
        tmp_iflim_res{:, 'check_max'} = tmp_iflim_res{:, 'mu_max'} .* (tmp_iflim_res{:, 'qty'} - tmp_iflim_res{:, 'max'});
        if isfield(mpc.if, 'name')
            tmp_iflim_res = [array2table(mpc.if.name, 'VariableNames', {'name'}), tmp_iflim_res];
        end
        iflim_res = [iflim_res; tmp_iflim_res];
    end
    % annual average results for each interface
    for i = 1:nc
        idx = iflim_res{:, 'if_num'} == i;
        tmp_res = iflim_res(idx, :);
        tmp_iflim_res = tmp_res;
        tmp_iflim_res(2:end, :) = [];
        tmp_iflim_res{:, 'hr_num'} = 0;
        tmp_iflim_res{:, 'min'} = min(tmp_res{:, 'min'});
        tmp_iflim_res{:, 'max'} = max(tmp_res{:, 'max'});
        tmp_iflim_res{:, 'qty'} = tmp_res{:, 'qty'}' * esc.hrs_map{:, 'probability'};
        tmp_iflim_res{:, 'mu_min'} = tmp_res{:, 'mu_min'}' * esc.hrs_map{:, 'probability'};
        tmp_iflim_res{:, 'mu_max'} = tmp_res{:, 'mu_max'}' * esc.hrs_map{:, 'probability'};
        tmp_iflim_res{:, 'mu'} = tmp_res{:, 'mu'}' * esc.hrs_map{:, 'probability'};
        tmp_iflim_res{:, 'flow_mwh'} = sum(tmp_res{:, 'flow_mwh'});
        tmp_iflim_res{:, 'merch_surplus'} = sum(tmp_res{:, 'merch_surplus'});
        tmp_iflim_res{:, 'check_min'} = 0; % Dont check annual average
        tmp_iflim_res{:, 'check_max'} = 0; % Dont check annual average

        iflim_res = [iflim_res; tmp_iflim_res];
    end

    res_constraint.iflim = iflim_res;
end



