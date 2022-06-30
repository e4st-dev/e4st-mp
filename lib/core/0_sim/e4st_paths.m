function eopt = e4st_paths(eopt)
%% E4ST_PATHS: Adds E4ST directories to the search path
% Sets default output folders if not pre-specified by user
% And creates output folders if not pre-existing

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Steven Witkin (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% FILE PATHS
% Get Main E4ST Directory
if ~isfield(eopt, 'root_dir') 
    cur_fn = fileparts(mfilename('fullpath'));
    eopt.root_dir = fullfile(cur_fn, '..', '..', '..');
end

% Add all subfolders of lib to the search path
if exist('e4st_solve', 'file') ~= 2 % Only need to do if not already on path
    lib_dir = fullfile(eopt.root_dir, 'lib');
    addpath(genpath(lib_dir))
end

if ~isfield(eopt, 'out_dir_res') % For outputs
    eopt.out_dir_res = fullfile(eopt.root_dir, 'Output', 'Results', eopt.prj_name);
end

if ~isfield(eopt, 'out_dir_raw') % A place to save raw results
    eopt.out_dir_raw = fullfile(eopt.root_dir, 'Output', 'RawRes', eopt.prj_name);
end

if ~isfield(eopt, 'out_dir_setup') % A place to save setup information
    eopt.out_dir_setup = fullfile(eopt.root_dir, 'Output', 'Setup', eopt.prj_name);
end

if ~isfield(eopt, 'out_dir_log') % For logs
    eopt.out_dir_log = fullfile(eopt.root_dir, 'Output', 'Logs', eopt.prj_name);
end

% Create directories if they do not already exist
for fol = {eopt.out_dir_log, eopt.out_dir_setup, eopt.out_dir_raw, eopt.out_dir_res}
    if exist(fol{:}, 'dir') ~= 7
        mkdir(fol{:})
    end
end

%% Check if MATPOWER is available
assert(exist('mpver', 'file') == 2, "MATPOWER must be on search path to run E4ST")

end