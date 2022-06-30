function mpc = updateGenCost(mpc, idx_gen, eopt)
% updateGenCost Update gencost ($/MWh) for specified generators
% - gencost = vom + (hr * fuel cost) - ptc
% - idx_gen specifies indices for generators to update

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Check generator index
if isempty(idx_gen)
    idx_gen = ones(size(mpc.gen), 1) == 1;
end

%% Do not update dispatchable load
idx_gen = idx_gen & ~strcmp(mpc.genfuel,'dl');

%% Update gencost: gencost = vom + (hr * fuel cost) - ptc
mpc.gencost(idx_gen, 5) = (mpc.gen_aux{idx_gen, 'FUEL_COST'} .* mpc.gen_aux{idx_gen,'HR'}) + ...
                           mpc.gen_aux{idx_gen, 'VOM'} - mpc.gen_aux{idx_gen, 'PTC'};

%% Display progress information
vfprintf(eopt.verbose, '(Gencost updated for %d generators)\n', sum(idx_gen));
