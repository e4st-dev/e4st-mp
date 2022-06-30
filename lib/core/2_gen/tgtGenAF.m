function mpc = tgtGenAF(mpc, esc, eopt)
%tgtGenAF: Target grouped generator availability factors

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Prepare Gen AF Target Data
info = esc.gen_tgt_af;
info = filterInfo(info, esc);

if isfield(info, 'joint')
    info_list = info.joint{:, :};
else
    info_list = (1:height(info.value));
end

if isfield(mpc, 'idx_update_af')
    idx_update_af = mpc.idx_update_af;
else
    idx_update_af = ones(size(mpc.gen, 1), 1) == 1;
end

%% Target Gen AF
define_constants;
probability = esc.hrs_map{:, 'probability'};
idx_online = mpc.gen(:, GEN_STATUS) == 1;

for i = unique(info_list)'
    idx_info = i == info_list;

    cur_info = filterStruct(info, idx_info);
    idx_gen = getInfoIdx(mpc, esc, cur_info, 'gen');
    idx_gen = idx_gen & idx_online & idx_update_af;

    if ~any(idx_gen)
        continue
    end

    tgt = unique(cur_info.value{:, 'tgt_af'});
    delta = 999;
    if cur_info.agg{:, :} == 0
        for i = find(idx_gen)
            while delta > .005
                avg_af = sum(mpc.availability_factor(i, :)*probability.*mpc.gen(i, PMAX)) / sum(mpc.gen(i, PMAX));
                mpc.availability_factor(i, :) = mpc.availability_factor(i, :) * (tgt / avg_af);
                mpc.availability_factor(mpc.availability_factor > 1) = 1;
                new_avg = sum(mpc.availability_factor(i, :)*probability.*mpc.gen(i, PMAX)) / sum(mpc.gen(i, PMAX));
                delta = abs(new_avg-avg_af);
            end
        end
    else
        while delta > .005
            avg_af = sum(mpc.availability_factor(idx_gen, :)*probability.*mpc.gen(idx_gen, PMAX)) / sum(mpc.gen(idx_gen, PMAX));
            mpc.availability_factor(idx_gen, :) = mpc.availability_factor(idx_gen, :) * (tgt / avg_af);
            mpc.availability_factor(mpc.availability_factor > 1) = 1;
            new_avg = sum(mpc.availability_factor(idx_gen, :)*probability.*mpc.gen(idx_gen, PMAX)) / sum(mpc.gen(idx_gen, PMAX));
            delta = abs(new_avg-avg_af);
        end
    end

end

%% Target Gen AF
fprintf('Targeting generator AFs\n')
idx_tgt_af = find(idx_update_af & mpc.gen_aux{:, 'TGT_AF'} > 0);
tgt_afs = mpc.gen_aux{:, 'TGT_AF'};
if any(idx_tgt_af)
    for i = idx_tgt_af'
        tgt = tgt_afs(i);
        delta = 999;
        while delta > .005
            avg_af = sum(mpc.availability_factor(i, :)*probability.*mpc.gen(i, PMAX)) / sum(mpc.gen(i, PMAX));
            mpc.availability_factor(i, :) = mpc.availability_factor(i, :) * (tgt / avg_af);
            mpc.availability_factor(mpc.availability_factor > 1) = 1;
            new_avg = sum(mpc.availability_factor(i, :)*probability.*mpc.gen(i, PMAX)) / sum(mpc.gen(i, PMAX));
            delta = abs(new_avg-avg_af);
        end
        vfprintf(eopt.verbose, 'Avg AF %s generator set to %.2f (target: %.2f)\n', mpc.gentype{i}, new_avg, tgt);
    end
end