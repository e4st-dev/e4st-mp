function [mpc, result] = RunE4ST(mpc, esc, pre_result, eopt)
%RunE4ST: Run a single-year E4ST simulation

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% ----------E4ST SIM OPTS--------------
eopt = e4st_opts(eopt);
eopt = e4st_paths(eopt);

%% START LOG
%create log directory
log_dir = fullfile(eopt.out_dir_log, eopt.case_name, eopt.grid_name);
if ~exist(log_dir, 'dir')
    mkdir(log_dir)
end
%convert eopt verbose to structure, if not already
if ~isstruct(eopt.verbose)
    verbose = eopt.verbose;
    eopt = rmfield(eopt, 'verbose');
    eopt.verbose.verbose = verbose;
end
%specify filename to which to print
if eopt.iter.status
   eopt.verbose.filename = fullfile(log_dir, [eopt.sim_name_yr_iter, '.txt']);
else
   eopt.verbose.filename = fullfile(log_dir, [eopt.sim_name_yr, '.txt']);
end


%% ----------SETUP E4ST SIMULATION --------------
tic % Time simulation
if eopt.load_setup.status
    fn = fullfile(eopt.out_dir_setup, eopt.case_name, [eopt.sim_name_yr, '.mat']);
    load(fn, 'mpc', 'esc', 'offer', 'contab', 'eopt')
    vfprintf(eopt.verbose, 'Loaded E4ST inputs from %s \n', fn)
else
    [mpc, esc, offer, contab, eopt] = e4st_setup(mpc, esc, eopt);
end
% Save E4ST Inputs
if eopt.save_setup.status
    out_folder = fullfile(eopt.out_dir_setup, eopt.case_name);
    if exist(out_folder, 'dir') ~= 7
        mkdir(out_folder)
    end
    save(fullfile(out_folder, [eopt.sim_name_yr, '.mat']), 'mpc', 'esc', 'offer', 'contab', 'eopt') %'-v7.3'
    vfprintf(eopt.verbose, 'Saved E4ST inputs to %s \n', out_folder)
end
setup_time = toc; tic

%% ----------CALL CORE E4ST CODE-----------------
res = e4st_core(mpc, offer, contab, eopt);
run_time = toc; tic

%% ----------POST-PROCESS E4ST RESULTS-----------
[mpc, result] = e4st_post(res, mpc, esc, offer, contab, pre_result, eopt);
post_time = toc;

%% Save Run Time Summary (NEEDS CLEANUP)
if eopt.iter.status
    rt = array2table({eopt.sim_name_yr_iter}, 'VariableNames', {'sim_name_yr'});
else
    rt = array2table({eopt.sim_name_yr}, 'VariableNames', {'sim_name_yr'});
end
%rt = array2table({eopt.sim_name_yr}, 'VariableNames',{'sim_name_yr'});
rt{:, 'setup_time'} = setup_time / 60;
rt{:, 'run_time'} = run_time / 60;
rt{:, 'postprocess_time'} = post_time / 60;
rt{:, 'total_time'} = (setup_time + run_time + post_time) / 60;
result.(eopt.sim_yr).user_result.slv_res.run_time = rt;
if eopt.iter.status
    result.(eopt.sim_yr_iter).user_result.slv_res.run_time = rt;
end

out_folder = fullfile(eopt.out_dir_res, eopt.case_name, eopt.grid_name);
out_file = fullfile(out_folder, [eopt.sim_name_yr, '_runtime.xlsx']);
writetable(rt, out_file, 'Sheet', 'runtime')

%% Save Sim Year Names
sim_yrs = [];
for field = fieldnames(result)'
    if strcmp(field, 'inputs')
        continue
    end
    sim_yrs = [sim_yrs; field];
end
out_folder = fullfile(eopt.out_dir_res, eopt.case_name, eopt.grid_name);
out_file = fullfile(out_folder, [eopt.case_grid, '_simyrs.xlsx']);
writetable(table(sim_yrs), out_file, 'Sheet', 'sim_yrs')

%% Aggregate Case (Grid-Year) Results
fprintf('\nAggregating Case (Grid-Year) Results\n')
all_res.(eopt.case_grid) = result;

all_res.([eopt.case_name, '_all']) = result_user_agg(all_res, eopt);
tmp = all_res.([eopt.case_name, '_all']);
clear all_res;
all_res.([eopt.case_name, '_all']) = tmp;

out_folder = fullfile(eopt.out_dir_res, eopt.case_name, eopt.grid_name, 'all');
if exist(out_folder, 'dir') ~= 7
    mkdir(out_folder)
end

out_file = fullfile(out_folder, [[eopt.case_name, '_all'], '.mat']);
save(out_file, 'all_res')

out_file = fullfile(out_folder, [[eopt.case_name, '_all'], '.xlsx']);
if isfile(out_file)
    delete(out_file)
end

subfields = fieldnames(all_res.([eopt.case_name, '_all']));
for subfield = subfields'
    if isa(all_res.([eopt.case_name, '_all']).(char(subfield)), 'double')
        writetable(array2table(all_res.([eopt.case_name, '_all']).(char(subfield))), out_file, 'Sheet', char(subfield))
    elseif isa(all_res.([eopt.case_name, '_all']).(char(subfield)), 'table')
        writetable(all_res.([eopt.case_name, '_all']).(char(subfield)), char(out_file), 'Sheet', char(subfield))
    end
end

%% Run Key Results
fprintf('\nWriting Key Results \n')
if isfield(eopt, 'key_areas')
    areas = eopt.key_areas;
else
    for area = eopt.rpt_areas
        areas.(char(area)) = [];
    end
end
key_res_long = result_user_agg_key(all_res.([eopt.case_name, '_all']), areas, eopt);

out_file = fullfile(out_folder, [[eopt.case_name, '_key_res_long'], '.mat']);
save(out_file, 'key_res_long')

out_file = fullfile(out_folder, [[eopt.case_name, '_key_res_long'], '.xlsx']);
if isfile(out_file)
    delete(out_file)
end
fields = fieldnames(key_res_long);
for field = fields'
    writetable(key_res_long.(char(field)), out_file, 'Sheet', char(field))
end

%% Report Run Time Summary
vfprintf(eopt.verbose, '\n---------------------------------------------------------\n');
if eopt.iter.status
    vfprintf(eopt.verbose, 'Simulation Finished: %s\n\n', eopt.sim_name_yr_iter)
else
    vfprintf(eopt.verbose, 'Simulation Finished: %s\n\n', eopt.sim_name_yr)
end
vfprintf(eopt.verbose, 'Run Time Summary:\n')
vfprintf(eopt.verbose, '\tSetup Time: %.3f seconds (%.1f minutes)\n', setup_time, setup_time/60)
vfprintf(eopt.verbose, '\tRun Time: %.3f seconds (%.1f minutes)\n', run_time, run_time/60)
vfprintf(eopt.verbose, '\tPost-Process Time: %.3f seconds (%.1f minutes)\n', post_time, post_time/60)
vfprintf(eopt.verbose, '\tTotal Time: %.3f seconds (%.1f minutes)\n', setup_time+run_time+post_time, (setup_time + run_time + post_time)/60)
vfprintf(eopt.verbose, '---------------------------------------------------------\n');

