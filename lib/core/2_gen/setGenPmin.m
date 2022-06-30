function [mpc, offer] = setGenPmin(mpc, esc, offer, eopt)
%setGenPmin

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Prepare Gen PMIN Data
info = esc.gen_pmin;
info = filterInfo(info, esc);

%% Set Gen PMN
define_constants;
for i = 1:height(info.value)
    idx_info = i == (1:height(info.value))';
    cur_info = filterStruct(info, idx_info);
    idx_gen = getInfoIdx(mpc, esc, cur_info, 'gen');

    if ~any(idx_gen)
        continue
    end

    %mpc.gen(idx_gen, PMIN) = mpc.gen(idx_gen, PMAX) * cur_info.value{:, :};
    offer(idx_gen, 13) = cur_info.value{:, :};

    vfprintf(eopt.verbose, 'Pmin set to %.2f for %d generators\n', cur_info.value{:, :}, sum(idx_gen));
end
