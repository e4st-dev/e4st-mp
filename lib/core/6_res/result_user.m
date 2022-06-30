function res = result_user(res, mpc, esc, offer, contab, eopt)
% result_user runs several result subfunctions

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

% Sim Info
res.user_result.sim_info = result_siminfo(esc, eopt);

% Detailed Results
res.user_result.bus_res = result_bus(res, mpc, esc, eopt);
res.user_result.gen_res = result_gen(res, mpc, offer, contab, esc, eopt);
res.user_result.pol_res.pol_info = res.user_result.gen_res.policy; %result_policy called inside result_gen
res.user_result.gen_res = rmfield(res.user_result.gen_res, 'policy');
res.user_result.branch_res = result_branch(res, mpc, esc);
res.user_result.constraints_res = result_constraints(res, mpc, esc);
if ~isfield(eopt, 'rerun_results') && strcmp(eopt.solver, 'gurobi') %if not reprocessing the results and used gurobi
    res.user_result.slv_res = result_solution(res);
else
    res.user_result.slv_res.status = [];
end

% Summarized Results by Area
if isfield(eopt, 'rpt_areas')
    areas = eopt.rpt_areas;
    res.user_result = result_byarea(res.user_result, esc, areas, eopt);
elseif isfield(esc, 'report_res')
    areas = esc.report_res.area{esc.report_res.status{:, :} == 1, :}';
    res.user_result = result_byarea(res.user_result, esc, areas, eopt);
end
