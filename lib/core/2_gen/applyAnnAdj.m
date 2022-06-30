function [mpc, offer] = applyAnnAdj(mpc, offer, esc, adj_type, eopt)
%applyAnnAdj: apply annual adjustments and set values that depend on year

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Pre-Process Data
info = esc.gen_ann_adj;
info = filterInfo(info, esc);

idx_info = strcmp(info.var_name{:, :}, adj_type);%get only columns with the appropraite adj_type
info = filterStruct(info, idx_info);

%info.value = info.value(:, ['Y' num2str(esc.year)]);
idx_year = find(strcmp(info.value.Properties.VariableNames, ['Y', num2str(esc.year)]));

idx_dl = strcmp(mpc.genfuel, 'dl');

%% Apply Changes
for i = 1:height(info.value)
    idx_info = i == (1:height(info.value))';
    cur_info = filterStruct(info, idx_info);
    idx_gen = getInfoIdx(mpc, esc, cur_info, 'gen');%index of relevant gens
    adj_var = cur_info.col_name{:, :};
    idx_gen = idx_gen & ~idx_dl;

    if ~any(idx_gen)
        continue
    end
    
    method = char(cur_info.method{:, :});

    ann_val = cur_info.value{:, idx_year};
    cum_val = sum(cur_info.value{:, (idx_year - esc.year_delta):idx_year});

    switch adj_type
        case 'gen_aux'
            adj_data = mpc.gen_aux{:, adj_var};
        case 'offer'
            adj_data = offer(:, str2double(adj_var));
        case 'gencost'
            adj_data = mpc.gencost(:, str2double(adj_var));
        case 'none'
            continue
        otherwise
            vfprintf(eopt.verbose, 'Adjustment dataset %s not recognized\n', adj_type)
    end

    switch method
        case 'scale'
            adj_data(idx_gen, 1) = adj_data(idx_gen, 1) * ann_val;
            vfprintf(eopt.verbose, 'Applied "%s" method for %d generators to %s with value %.2f\n', method, sum(idx_gen), char(adj_var), ann_val)
        case 'set'
            adj_data(idx_gen, 1) = ann_val;
            vfprintf(eopt.verbose, 'Applied "%s" method for %d generators to %s with value %.2f\n', method, sum(idx_gen), char(adj_var), ann_val)
        case 'add'
            adj_data(idx_gen, 1) = adj_data(idx_gen, 1) + ann_val;
            vfprintf(eopt.verbose, 'Applied "%s" method for %d generators to %s with value %.2f\n', method, sum(idx_gen), char(adj_var), cum_val)
        case 'add_cum'
            adj_data(idx_gen, 1) = adj_data(idx_gen, 1) + cum_val;
            vfprintf(eopt.verbose, 'Applied "%s" method for %d generators to %s with value %.2f\n', method, sum(idx_gen), char(adj_var), cum_val)
        otherwise
            vfprintf(eopt.verbose, 'Method "%s" not recognized\n', method)
    end

    %% Apply Changes
    switch adj_type
        case 'gen_aux'
            mpc.gen_aux{:, adj_var} = adj_data;
            mpc = updateGenCost(mpc, idx_gen, eopt);
            offer = updateOfferPrc(mpc, offer, idx_gen, eopt);
        case 'offer'
            offer(:, str2double(adj_var)) = adj_data;
        case 'gencost'
            mpc.gencost(:, str2double(adj_var)) = adj_data;
        otherwise
            vfprintf(eopt.verbose, 'Adjustment dataset %s not recognized\n', adj_type)
    end

end
