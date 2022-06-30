function key_res_long = result_user_agg_key(all_res, areas, opt)
% result_user_agg_key collects key results
% if opt.key_res_long = 'T', pivots results into a long format

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%%
opt.key_res = 'T';
opt.key_res_long = 'T';

%% Areas
all_areas = {'grid', 'interconnection', 'nation', 'state'}; % Keep consistent with result_byarea.m
for area = all_areas
    if ~isfield(areas, area)
        areas.(char(area)) = [];
    end
end

%% Key Results
if strcmp(opt.key_res, 'T')
    %out_dir = fullfile(root_dir, ['prj_' prj], sim_set);
    %out_folder = fullfile(out_dir, 'all');
    %out_file = fullfile(out_folder, [sim_set '_key_res']);
    key_res = [];
    %fields = fieldnames(all_res.([sim_set '_all']))';
    fields = {'bus', 'gen_on_gentype', 'pol_info', 'iflim', 'branch'};

    for field = fields
        field = char(field);
        fprintf('Processing result %s...', field)
        if ~isfield(all_res, field) %.([sim_set '_all'])
            fprintf('Does not exist\n')
            continue
        end
        data = all_res.(field); %.([sim_set '_all'])
        if ~any(strcmp(field, {'pol_info', 'iflim'}))
            idx = zeros(height(data), 1);
            for area = fieldnames(areas)'
                if isempty(areas.(char(area))) % all subareas
                    idx = idx | strcmp(data{:, 'area'}, area);
                else % subareas
                    for subarea = areas.(char(area))
                        if any(strcmp(field, {'branch'}))
                            idx = idx | (strcmp(data{:, 'area'}, area) & strcmp(data{:, 'tarea'}, subarea));
                            idx = idx | (strcmp(data{:, 'area'}, area) & strcmp(data{:, 'farea'}, subarea));
                        else
                            idx = idx | (strcmp(data{:, 'area'}, area) & strcmp(data{:, 'subarea'}, subarea));
                        end
                    end
                end
            end
            data = data(idx, :);
        end
        vars = data.Properties.VariableNames;
        data{:, 'casename'} = {opt.case_name};
        data = data(:, [{'casename'}, vars]);

        if any(strcmp(field, {'branch_ac', 'branch_dc'})) && ~isfield(key_res, 'branch')
            idx_var = find(strcmp(data.Properties.VariableNames, 'tarea'));
            if strcmp(field, 'branch_ac')
                data{:, 'br_type'} = {'ac'};
            elseif strcmp(field, 'branch_dc')
                data{:, 'br_type'} = {'dc'};
            end
            data = [data(:, 1:idx_var), data(:, 'br_type'), data(:, (idx_var + 1):(end -1))];
            field = 'branch';
            key_res.(char(field)) = data;
        elseif any(strcmp(field, {'branch_ac', 'branch_dc'})) && isfield(key_res, 'branch')
            idx_var = find(strcmp(data.Properties.VariableNames, 'tarea'));
            if strcmp(field, 'branch_ac')
                data{:, 'br_type'} = {'ac'};
            elseif strcmp(field, 'branch_dc')
                data{:, 'br_type'} = {'dc'};
            end
            data = [data(:, 1:idx_var), data(:, 'br_type'), data(:, (idx_var + 1):(end -1))];
            field = 'branch';
            key_res.(char(field)) = [key_res.(char(field)); data];
        else
            key_res.(char(field)) = data;
        end

        %writetable(data, [out_file '.xlsx'], 'Sheet', field)

        fprintf('Finished\n')
    end
end

%% KEY RES LONG
if strcmp(opt.key_res_long, 'T')
    %out_dir = fullfile(root_dir, ['prj_' prj], sim_set);
    %out_folder = fullfile(out_dir, 'all');
    %out_file = fullfile(out_folder, [sim_set '_key_res_long']);
    fields = {'bus', 'gen_on_gentype', 'pol_info', 'branch', 'iflim'};
    %for field = fieldnames(key_res)'
    for field = fields
        fprintf('Processing result %s...', char(field))
        if isfield(key_res, char(field))
            data = key_res.(char(field));
        else
            fprintf('Does not exist\n')
            continue
        end
        switch char(field)
            case 'bus'
                key_res_long.(char(field)) = [];
                idx_desc = find(strcmp(data.Properties.VariableNames, 'GroupCount'));
                idx_lat = find(strcmp(data.Properties.VariableNames, 'latitude'));
                desc = data(:, [1:idx_desc, idx_lat, (idx_lat + 1)]);
                %vals = data(:,(idx_desc+1):end);
                data(:, [1:idx_desc, idx_lat, (idx_lat + 1)]) = [];
            case 'gen_on_gentype'
                key_res_long.(char(field)) = [];
                idx_desc = find(strcmp(data.Properties.VariableNames, 'GroupCount'));
                idx_lat = find(strcmp(data.Properties.VariableNames, 'latitude'));
                desc = data(:, [1:idx_desc, idx_lat, (idx_lat + 1)]);
                data(:, [1:idx_desc, idx_lat, (idx_lat + 1)]) = [];
            case 'pol_info'
                key_res_long.(char(field)) = [];
                idx_desc = find(strcmp(data.Properties.VariableNames, 'pol_name'));
                desc = data(:, 1:idx_desc);
                data(:, 1:idx_desc) = [];
            case 'branch'
                key_res_long.(char(field)) = [];
                idx_desc = find(strcmp(data.Properties.VariableNames, 'GroupCount'));
                idx_hr = find(strcmp(data.Properties.VariableNames, 'C1'));
                desc = data(:, 1:idx_desc);
                data(:, [(1:idx_desc), (idx_hr:end)]) = [];
            case 'iflim'
                key_res_long.(char(field)) = data;
            case 'branch_ac'
            case 'branch_dc'
        end
        if ~any(strcmp(field, {'iflim'}))
            vars = data.Properties.VariableNames;
            for var = vars
                desc{:, 'varname'} = var;
                desc{:, 'value'} = data{:, var};
                key_res_long.(char(field)) = [key_res_long.(char(field)); desc];
            end
        end

        %writetable(key_res_long.(char(field)),[out_file '.xlsx'],'Sheet',char(field))

        fprintf('Finished\n')

    end
end
