function [mpc, offer] = append_e4st_gen(mpc, offer, new_gen)
% append_e4st_gen appends mpc generators with full info onto the main mpc and offer

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

fields = {'gen', 'gencost', 'genfuel', 'gentype', ...
    'newgen', 'gen_aux', 'gen_desc', 'gen_map'};
for field = fields
    mpc.(char(field)) = [mpc.(char(field)); new_gen.(char(field))];
end
offer = [offer; new_gen.offer];

%add storage if it exists
if isfield(mpc, 'short_term_storage') && isfield(new_gen, 'short_term_storage')
    mpc.short_term_storage = [mpc.short_term_storage; new_gen.short_term_storage];
elseif isfield(mpc, 'short_term_storage')
    mpc.short_term_storage = [mpc.short_term_storage; zeros(size(new_gen.gen,1), 3)];
else
    error('Field mpc.short_term_storage not found')
end

end

