function eopt = e4st_opts(eopt)
%% E4ST_OPTS: Sets default e4st options if not pre-specified by user
%  eopt: e4st options struct

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% MAIN OPTIONS
if ~isfield(eopt, 'verbose') % Print Progress in Consule
    eopt.verbose = 1;
end
if ~isfield(eopt, 'default') % "Default" Simulation (for PRL)
    eopt.default = 'F';
end
if ~isfield(eopt, 'prl') % hourly price responsive load
    eopt.prl = 'F';
end
if ~isfield(eopt, 'iter') % simulation involves iterations
    eopt.iter.status = 0;
    eopt.iter.name = '';
end
if ~isfield(eopt, 'update_gens') % apply sim updates to generators
    eopt.update_gens = 'T';
end

%% Solving options
if ~isfield(eopt, 'threads')
    eopt.threads = 0;
end
if ~isfield(eopt, 'solver') % 'linprog' or 'gurobi'
    eopt.solver = 'gurobi';
end
if ~isfield(eopt, 'NumericFocus') % Gurobi Numeric Focus option
    eopt.NumericFocus = 0;
end
if ~isfield(eopt, 'BarHomogeneous') % Gurobi BarHomogeneous option
    eopt.BarHomogeneous = -1;
end

%% Generator Options
if ~isfield(eopt, 'ccs') % Existence of CCS
    eopt.ccs = 'F';
end
if ~isfield(eopt, 'dac') % Existence of DAC
    eopt.dac = 'F';
end
if ~isfield(eopt, 'ccs_retro') % CCS Retrofits
    eopt.ccs_retro = 'F';
end
if ~isfield(eopt, 'hydrogen_retro') % Hydrogen Retrofits on NG
    eopt.hydrogen_retro = 'F';
end
if ~isfield(eopt, 'storage') % Battery Storage
    eopt.storage = 'F';
end

%% Input and Output Options
if ~isfield(eopt, 'save_raw_all')
    eopt.save_raw_all = 'F';
end
if ~isfield(eopt, 'save_hourly')
    eopt.save_hourly = 'F';
end
if ~isfield(eopt, 'save_setup') %saves inputs into e4st_solve for future reference
    eopt.save_setup.status = 1;
end
if ~isfield(eopt, 'load_setup')
    eopt.load_setup.status = 0;
end
if ~isfield(eopt, 'sum_res') % run detailed summaries of key results (default: T)
    eopt.sum_res = 'T';
end
if ~isfield(eopt, 'res2xlsx') % save key results as .xslx (default: T)
    eopt.res2xlsx = 'T';
end
if ~isfield(eopt, 'res2mat') % save key results as .mat (default: T)
    eopt.res2mat = 'T';
end
if ~isfield(eopt, 'raw2mat') % save extensive raw results as .mat (default: F)
    eopt.raw2mat = 'F';
end

%% Constants
if ~isfield(eopt, 'CH4_GWP_20y')
    eopt.CH4_GWP_20y = 86;
end

%% NAME CONVENTIONS
assert(isfield(eopt, 'case_name'), "At minimum, eopt.case_name must be defined") 

if ~isfield(eopt, 'prj_name') % Project name, for subfolders of results, etc
    eopt.prj_name = '';
end

if ~isfield(eopt, 'grid_name') % Grid name
    eopt.grid_name = '';
end

if ~isfield(eopt, 'sim_yr') % Actual year kept in esc.year, this is for naming purposes
    eopt.sim_yr = '';
end


% standard way of defining 'case_grid' based on 'case_name' and 'grid_name'
if isfield(eopt, 'case_name') && isfield(eopt, 'grid_name')
    eopt.case_grid = [eopt.case_name, '_', eopt.grid_name];
end

%set simulation name in a standard way
if strcmp(eopt.default, 'T') && isfield(eopt, 'case_grid') % Default simulation: First simulation year, no PRL
    eopt.sim_name = ['default_', eopt.case_grid];
elseif isfield(eopt, 'case_grid')
    eopt.sim_name = ['result_', eopt.case_grid];
end

%combined simulation year and name
eopt.sim_name_yr = [eopt.sim_name, '_', eopt.sim_yr];

if eopt.iter.status % if iterations turned on
    eopt.sim_yr_iter = [eopt.sim_yr, '_i', eopt.iter.name];
    eopt.sim_name_yr_iter = [eopt.sim_name, '_', eopt.sim_yr_iter];
end
