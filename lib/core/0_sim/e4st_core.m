function res = e4st_core(mpc, offer, contab, eopt)
% e4st_core Set options and run e4st_solve

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Gurobi OPTIONS
% https://www.gams.com/latest/docs/S_GUROBI.html
% TERMINATION OPTIONS
%bariterlimit	Limits the number of barrier iterations performed	infinity
%cutoff	Sets a target objective value	0
%iterationlimit	Limits the number of simplex iterations performed	infinity
%nodelimit	Limits the number of MIP nodes explored	maxdouble
%solutionlimit	Limits the number of feasible solutions found	maxint
%timelimit	Limits the total time expended in seconds	GAMS reslim

% TOLERANCE OPTIONS
%barconvtol	Controls barrier termination	1e-8
%barqcpconvtol	Convergence tolerance for the barrier algorithm when solving a QCP	1e-6
%feasibilitytol	Primal feasibility tolerance	1e-6
%intfeastol	Integer feasibility tolerance	1e-5
%markowitztol	Threshold pivoting tolerance	0.0078125
%mipgap	Relative MIP optimality gap	GAMS optcr
%mipgapabs	Absolute MIP optimality gap	GAMS optca
%optimalitytol	Dual feasibility tolerance	1e-6
%psdtol	limit on the amount of diagonal perturbation	1e-6

% SIMPLEX OPTIONS
%normadjust	Pricing norm variants	-1
%objscale	Objective coefficients scaling	0
%perturbvalue	Magnitude of simplex perturbation when required	0.0002
%quad	Quad precision computation in simplex	-1
%scaleflag	Enables or disables model scaling	1
%sifting	Sifting within dual simplex	-1
%siftmethod	LP method used to solve sifting sub-problems	-1
%simplexpricing	Determines variable pricing strategy	-1

% BARRIER OPTIONS
%barcorrectors	Limits the number of central corrections performed in each barrier iteration	-1
%barhomogeneous	Homogeneous barrier algorithm selection	-1
%barorder	Chooses the barrier sparse matrix fill-reducing algorithm	-1
%crossover	Determines the crossover strategy used to transform the barrier solution into a basic solution	-1
%crossoverbasis	Determines the initial basis construction strategy for crossover	0
%qcpdual	Determines whether dual variable values are computed for QCP models	1

% Gurobi Output
% the progress of the barrier algorithm in iterations with the primal and dual objective values,
% the magnitude of the primal and dual infeasibilites,
% and the magnitude of the complementarity violation.

%% MATPOWER OPTIONS
threads = eopt.threads;

mpopt = mpoption('model', 'DC');
mpopt = mpoption(mpopt, 'out.all', 0);
mpopt = mpoption(mpopt, 'verbose', 3);
switch eopt.solver
    case 'gurobi'
        mpopt = mpoption(mpopt, 'gurobi.threads', threads);
        mpopt = mpoption(mpopt, 'gurobi.method', 2); % Turn on barrier method
        mpopt = mpoption(mpopt, 'gurobi.opts.BarIterLimit', 1000); % Max number of barrier iterations
        mpopt = mpoption(mpopt, 'gurobi.opts.Crossover', 0); % Turn off crossover
        mpopt = mpoption(mpopt, 'gurobi.opts.FeasibilityTol', 1e-02); % Primal Tolerance (Column 3)
        mpopt = mpoption(mpopt, 'gurobi.opts.OptimalityTol', 1e-06); % Dual Tolerance (Column 4)
        mpopt = mpoption(mpopt, 'gurobi.opts.BarConvTol', 1e-06); % Complementarity Tolerance (Column 5)
        mpopt = mpoption(mpopt, 'gurobi.opts.LogFile', ...
            fullfile(eopt.out_dir_log, eopt.case_name, eopt.grid_name, 'gurobi.txt'));
        % mpopt = mpoption(mpopt, 'gurobi.opts.FeasibilityTol', 1e-02); % Primal Tolerance (Column 3)
        % mpopt = mpoption(mpopt, 'gurobi.opts.OptimalityTol',  1e-02);  % Dual Tolerance (Column 4)
        % mpopt = mpoption(mpopt, 'gurobi.opts.BarConvTol', .9);     % Complementarity Tolerance (Column 5)
        mpopt = mpoption(mpopt, 'gurobi.opts.NumericFocus', eopt.NumericFocus); % tradeoff between accuracy and solving (0 is default, 2 is focused)
        mpopt = mpoption(mpopt, 'gurobi.opts.BarHomogeneous', eopt.BarHomogeneous); % homogenous barrier method for detecting infeasabilities (-1 is default, 1 is on)       
        %mpopt = gurobi_user_options(mpopt, eopt);
        %if isfield(eopt, 'grb_opt')
        %end
    case 'linprog'
        mpoption(mpopt, 'opf.dc.solver', 'linprog_ds');
end

%% Assertions
% Check for missing values

assert(~any(isnan(contab), 'all'), 'Missing value in contab')
assert(~any(isnan(offer), 'all'), 'Missing value in offer')
assert(~any(isnan(mpc.bus), 'all'), 'Missing value in mpc.bus')
assert(~any(isnan(mpc.branch), 'all'), 'Missing value in mpc.branch')
assert(~any(isnan(mpc.dcline), 'all'), 'Missing value in mpc.dcline')
assert(~any(isnan(mpc.gen), 'all'), 'Missing value in mpc.gen')
assert(~any(isnan(mpc.gencost), 'all'), 'Mising value in mpc.gencost')
assert(~any(isnan(mpc.total_output.map), 'all'), 'Missing value in mpc.total_output.map')
assert(~any(isnan(mpc.total_output.max), 'all'), 'Missing value in mpc.total_output.max')
assert(~any(isnan(mpc.total_output.min), 'all'), 'Missing value in mpc.total_output.min')
assert(~any(isnan(mpc.total_output.coeff), 'all'), 'Missing value in mpc.total_output.coeff')
%assert(~any(isnan(mpc.caplim.max), 'all'), 'Missing value in mpc.caplim.max')
%assert(~any(isnan(mpc.caplim.min), 'all'), 'Missing value in mpc.caplim.min')
%assert(~any(isnan(mpc.caplim.map), 'all'), 'Missing value in mpc.caplim.map')
assert(~any(isnan(mpc.short_term_storage), 'all'), 'Missing value in mpc.short_term_storage')
assert(~any(isnan(mpc.availability_factor), 'all'), 'Missing value in mpc.availability_factor')
assert(~any(isnan(mpc.hourly_shape), 'all'), 'Missing value in mpc.hourly_shape')
assert(~any(isnan(mpc.hourly_load), 'all'), 'Missing value in mpc.hourly_load')

%% SOLVE E4ST
res = e4st_solve(mpc, offer, contab, mpopt);
