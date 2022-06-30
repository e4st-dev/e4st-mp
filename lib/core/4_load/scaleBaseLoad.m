function mpc = scaleBaseLoad(mpc, esc, eopt)
% scaleBaseLoad: Scale base hour loads in each load area

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Paul Picciano (Resources for the Future) and Biao Mao (Rensselaer Polytechnic Institute)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Prepare Load Growth Data
info = esc.load_growth;
info = filterInfo(info, esc);

info.value = info.value(:, ['Y', num2str(esc.year)]);

%% Scale base hour load
opt = struct('pq', 'P');
for i = 1:height(info.value)
    idx_info = i == (1:height(info.value))';
    cur_info = filterStruct(info, idx_info);
    idx_bus = getInfoIdx(mpc, esc, cur_info, 'bus');

    if ~any(idx_bus)
        continue
    end

    load_growth = cur_info.value{:, :}.^esc.year_delta;
    [mpc.bus, mpc.gen] = scale_load(load_growth, mpc.bus, mpc.gen, idx_bus, opt);
    vfprintf(eopt.verbose, 'Base loads scaled by %.2f in %d buses\n', load_growth, sum(idx_bus));
end
