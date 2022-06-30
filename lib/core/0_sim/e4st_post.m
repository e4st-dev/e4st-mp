function [mpc, result] = e4st_post(res, mpc, esc, offer, contab, pre_result, eopt)
% e4st_post processes and saves E4ST results

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).


%% Retire endogenously shutdown generators
if esc.year >= esc.first_ret_yr
    mpc = retireGen(mpc, offer, esc, res, eopt);
end

%%
if isfield(res.opf_results.raw.output, 'status') 
    res.slv_status = res.opf_results.raw.output.status;
else
    res.slv_status = 'na';
end

%% Create out folder
out_folder = fullfile(eopt.out_dir_res, eopt.case_name, eopt.grid_name);
if exist(out_folder, 'dir') ~= 7
    mkdir(out_folder)
end
%% Save mat results in case there is error in user results
% Instead of saving after, in the res2mat block below.
if strcmp(eopt.res2mat, 'T')
    vfprintf(eopt.verbose, 'Save Non-User Results to .MAT File\n')
    % mpc, offer, contab, esc, and eopt should not change from call to
    % result_user, so they are not saved again below
    out_file = fullfile(out_folder, [eopt.sim_name_yr, '_mp.mat']);
    save(out_file, 'mpc', 'offer', 'contab');
    out_file = fullfile(out_folder, [eopt.sim_name_yr, '_e4st.mat']);
    save(out_file, 'esc', 'eopt');
    % Saving doesn't work with full opf results.
    opf_results = res.opf_results;
    res.opf_results = [];
    res.opf_results.f = opf_results.f;
    out_file = fullfile(out_folder, [eopt.sim_name_yr, '_res.mat']);
    save(out_file, 'res');
    res.opf_results = opf_results;
    clear opf_results
end


%% Generate and Export User Results
if strcmp(eopt.sum_res, 'T')
    vfprintf(eopt.verbose, 'Generating User Results\n')
    res = result_user(res, mpc, esc, offer, contab, eopt);
    
    if strcmp(eopt.res2mat, 'T')
        vfprintf(eopt.verbose, 'Save User Results to .MAT File\n')
        
        %user_result = res.user_result;

        if strcmp(eopt.raw2mat, 'T')
           user_result.e4st.year = esc.year;
           user_result.e4st.esc = esc;
           user_result.e4st.mpc = mpc;
           user_result.e4st.offer = offer;
           user_result.e4st.contab = contab;
           user_result.e4st.eopt = eopt;
        end

        fVal = res.opf_results.f;
        res.opf_results = [];
        res.opf_results.f = fVal;
        out_file = fullfile(out_folder, [eopt.sim_name_yr, '_res.mat']);
        save(out_file, 'res');
        ur = res.user_result;
        out_file = fullfile(out_folder, [eopt.sim_name_yr, '_ur.mat']);
        save(out_file, 'ur');
        

        if eopt.iter.status
            out_file = fullfile(out_folder, [eopt.sim_name_yr_iter, '_res.mat']);
            save(out_file, 'res');
            out_file = fullfile(out_folder, [eopt.sim_name_yr_iter, '_mp.mat']);
            save(out_file, 'mpc', 'offer', 'contab');
            out_file = fullfile(out_folder, [eopt.sim_name_yr_iter, '_e4st.mat']);
            save(out_file, 'esc', 'eopt');
        end
    end
    
    if strcmp(eopt.res2xlsx, 'T')
        vfprintf(eopt.verbose, 'Save User Results to .XLSX File\n')
        fields = fieldnames(res.user_result); % {'bus_res', 'gen_res', 'branch', 'area_res', 'constraints_res'}
        for field = fields'
            subfields = fieldnames(res.user_result.(char(field)));
            for subfield = subfields'
                if (strcmp(field, 'bus_res') || strcmp(field, 'gen_res')) && ~strcmp(subfield, 'annual') %only write annual results to xlsx
                    continue
                end
                tmp_res = res.user_result.(char(field)).(char(subfield));
                if isempty(tmp_res)
                    continue
                end

                try
                    if eopt.iter.status
                        out_file = fullfile(out_folder, [eopt.sim_name_yr_iter, '_', char(field), '.xlsx']);
                    else
                        out_file = fullfile(out_folder, [eopt.sim_name_yr, '_', char(field), '.xlsx']);
                    end
                    if isa(tmp_res, 'double')
                        writetable(array2table(tmp_res), out_file, 'Sheet', char(subfield))
                    elseif isa(tmp_res, 'table')
                        writetable(tmp_res, out_file, 'Sheet', char(subfield))
                    end
                catch
                    if eopt.iter.status
                        out_file = fullfile(out_folder, [eopt.sim_name_yr_iter, '_', char(field), '_', datestr(now, 'HH_MM_SS_FFF'), '.xlsx']);
                    else
                        out_file = fullfile(out_folder, [eopt.sim_name_yr, '_', char(field), '_', datestr(now, 'HH_MM_SS_FFF'), '.xlsx']);
                    end
                    if isa(tmp_res, 'double')
                        writetable(array2table(tmp_res), out_file, 'Sheet', char(subfield))
                    elseif isa(tmp_res, 'table')
                        writetable(tmp_res, out_file, 'Sheet', char(subfield))
                    end
                end
            end
        end
    else
        vfprintf(eopt.verbose, 'Save Area Results to .XLSX File\n')
        fields = {'area_res'};
        for field = fields'
            subfields = fieldnames(res.user_result.(char(field)));
            for subfield = subfields'
                tmp_res = res.user_result.(char(field)).(char(subfield));
                if isempty(tmp_res)
                    continue
                end

                try
                    if eopt.iter.status
                        out_file = fullfile(out_folder, [eopt.sim_name_yr_iter, '_', char(field), '.xlsx']);
                    else
                        out_file = fullfile(out_folder, [eopt.sim_name_yr, '_', char(field), '.xlsx']);
                    end
                    if isa(tmp_res, 'double')
                        writetable(array2table(tmp_res), out_file, 'Sheet', char(subfield))
                    elseif isa(tmp_res, 'table')
                        writetable(tmp_res, out_file, 'Sheet', char(subfield))
                    end
                catch
                    if eopt.iter.status
                        out_file = fullfile(out_folder, [eopt.sim_name_yr_iter, '_', char(field), '_', datestr(now, 'HH_MM_SS_FFF'), '.xlsx']);
                    else
                        out_file = fullfile(out_folder, [eopt.sim_name_yr, '_', char(field), '_', datestr(now, 'HH_MM_SS_FFF'), '.xlsx']);
                    end
                    if isa(tmp_res, 'double')
                        writetable(array2table(tmp_res), out_file, 'Sheet', char(subfield))
                    elseif isa(tmp_res, 'table')
                        writetable(tmp_res, out_file, 'Sheet', char(subfield))
                    end
                end
            end
        end
    end



else
    res.user_result = [];
end

%% save hourly results if desired

if strcmp(eopt.save_hourly, 'T') && eopt.iter.status
    %set up hourly table .. similar to what is done in result_byarea
    hrs = esc.hrs_map{:, 'rep_hr'}';
    gen_hrs = [];
    lmp_hrs = [];
    for hr = hrs
        gen_hrs = [gen_hrs, {strjoin([{'GEN_'}, hr], '')}];
        lmp_hrs = [lmp_hrs, {strjoin([{'LMP_'}, hr], '')}];
    end
    ann_res = res.user_result.gen_res.annual;
    gen_hourly = array2table(res.user_result.gen_res.gen_hourly, 'VariableNames', gen_hrs);
    lmp_hourly = array2table(res.user_result.gen_res.lmp_hourly, 'VariableNames', lmp_hrs);
    hourly = [ann_res(:, {'genfuel', 'gentype', 'capacity_active_mw', 'generation_mwh', 'capacity_factor', 'lmp_gen_mwh', 'capacity_active_mw', 'capacity_start_mw'}), gen_hourly, lmp_hourly];
                
    %save results
    out_folder = fullfile(eopt.out_dir_res, eopt.case_name, eopt.grid_name);
    %save matlab file
    out_file = fullfile(out_folder, [eopt.sim_name_yr_iter, '_hourly_genlmp.mat']);
    save(out_file, 'hourly'); 
    %save excel file
    out_file = fullfile(out_folder, [eopt.sim_name_yr_iter, '_hourly_genlmp.xlsx']);
    writetable(hourly, out_file, 'Sheet', 'hourly_data')
end


%% Remove undesired results to reduce results size
fVal = res.opf_results.f;
res.opf_results = [];
res.opf_results.f = fVal;

%% Save Default Load and LMPs
if strcmp(eopt.default, 'T')
    mpc.default_price = res.user_result.bus_res.lmp_hourly;
    mpc.default_load = res.user_result.bus_res.load_hourly;
end

%% Save Results
result = pre_result;
result.(eopt.sim_yr) = res;
%result.(eopt.sim_yr).user_result = res.user_result;
result.(eopt.sim_yr).year = esc.year;
result.(eopt.sim_yr).e4st.esc = esc;
result.(eopt.sim_yr).e4st.mpc = mpc;
%result.(eopt.sim_yr).e4st.offer = offer;
%result.(eopt.sim_yr).e4st.contab = contab;
%result.(eopt.sim_yr).e4st.eopt = eopt;

if eopt.iter.status
    result.(eopt.sim_yr_iter) = result.(eopt.sim_yr);
end

if strcmp(eopt.save_raw_all, 'T')
    save([eopt.out_dir_raw, eopt.case_grid, '.mat'], 'result', '-v7.3');
end


