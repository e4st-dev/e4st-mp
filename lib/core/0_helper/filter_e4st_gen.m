function [mpc_subset, offer_subset] = filter_e4st_gen(mpc, offer, idx_gen)
% filter_e4st_gen returns a subset of mpc generator info based on idx_gen, 
% which is a logical vector that indexes rows of the mpc

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

fields = {'gen', 'gencost', 'genfuel', 'gentype', ...
    'newgen', 'gen_aux', 'gen_desc', 'gen_map'};
for field = fields
    data = mpc.(char(field));
    mpc_subset.(char(field)) = data(idx_gen, :);
end
offer_subset = offer(idx_gen, :);
end