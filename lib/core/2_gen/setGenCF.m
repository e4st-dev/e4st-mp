function mpc = setGenCF(mpc, esc, eopt)
% setGenCF: set generator capacity factors

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Prepare Gen CF Data
info = esc.gen_cf;
info = filterInfo(info, esc);

if isfield(info, 'joint')
    info_list = info.joint{:, :};
else
    info_list = (1:height(info.value));
end

%% Set Gen CF
define_constants;
[map, min, max, name] = deal([]);
idx_online = mpc.gen(:, GEN_STATUS) == 1;
idx_tgtcf = mpc.gen_aux{:, 'TGT_CF'};

for i = unique(info_list)'
    idx_info = i == info_list;

    cur_info = filterStruct(info, idx_info);
    idx_gen = getInfoIdx(mpc, esc, cur_info, 'gen');
    idx_gen = idx_gen & idx_online;

    if ~any(idx_gen)
        continue
    end

    min_val = unique(cur_info.value{:, 'min'});
    max_val = unique(cur_info.value{:, 'max'});

    if min_val == 0
        min_val = -Inf;
    end

    if max_val == 1
        max_val = Inf;
    end

    if cur_info.agg{:, :} == 0
        % map = idx_gen .* eye(size(mpc.gen, 1))
        % map = map(idx_gen,:);
        idx_gen = idx_gen & ~idx_tgtcf;
        mpc.gen_aux{idx_gen, 'TGT_CF'} = max_val;

        %tmp = zeros(1, size(mpc.gen, 1));
        %gens = find(idx_gen);
        %for j = 1:length(gens)
        %    tmp_map = tmp;
        %    tmp_map(1, gens(j)) = 1;
        %    map = [map; tmp_map];
        %    name = [name; {['CF_' num2str(i) '_' num2str(j)]}];
        %end
        %min = [min; mpc.gen(idx_gen, PG) * min_val];
        %max = [max; mpc.gen(idx_gen, PG) * max_val];
    elseif cur_info.agg{:, :} == 1
        map = [map; idx_gen'];
        name = [name; {['CF_Agg_', num2str(i)]}];
        min = [min; sum(mpc.gen(idx_gen, PG)) * min_val];
        max = [max; sum(mpc.gen(idx_gen, PG)) * max_val];
    end
    vfprintf(eopt.verbose, 'CF set to min of %.2f and max of %.2f for %d generators\n', cur_info.value{1, 1}, cur_info.value{1, 2}, sum(idx_gen));
end

%%
fprintf('Targeting generator CFs\n')
% Make sure not infeasible
hours = esc.hrs_map{:, 'hours'};
max_cf_limit = mpc.availability_factor * hours / sum(hours);
idx = mpc.gen_aux{:, 'TGT_CF'} > 0 & mpc.gen_aux{:, 'TGT_CF'} > max_cf_limit;
mpc.gen_aux{idx, 'TGT_CF'} = max_cf_limit(idx);

% Making a constraint for an offline (or PG = 0) gen will result in NaNs
idx_tgt_cf = mpc.gen_aux{:, 'TGT_CF'} > 0 & idx_online;
tgt_cfs = mpc.gen_aux{:, 'TGT_CF'};
gen_ids = mpc.gen_map{:, 'gen_id'};
%tmp = zeros(1, size(mpc.gen, 1));
%tic
if any(idx_tgt_cf)
    %     tmp_map = idx_tgt_cf .* eye(mpc.ng, mpc.ng);
    %     tmp_map = tmp_map(idx_tgt_cf, :);

    tmp_map = zeros(sum(idx_tgt_cf), mpc.ng);
    curr_row = 1;


    for i = find(idx_tgt_cf)'
        tmp_map(curr_row, i) = 1;
        tgt = tgt_cfs(i);
        id = gen_ids(i);
        %tmp_map = tmp;
        %tmp_map(1, i) = 1;
        %map = [map; tmp_map];
        name = [name; {['CF_Unagg_', num2str(id)]}];
        min = [min; mpc.gen(i, PG) * (tgt * -Inf)]; %
        max = [max; mpc.gen(i, PG) * (tgt * 1)];
        curr_row = curr_row + 1;
        vfprintf(eopt.verbose, 'Max CF for %s generator set to %.2f\n', mpc.gentype{i}, tgt);
    end
    map = [map; tmp_map];
end
%toc

%% Add TOC
coeff = ones(size(mpc.gen, 1), 1);
type = ones(length(max), 1);
mpc = addTOC(mpc, map, min, max, coeff, type, name);
