function mpc = applyPTC(mpc, credit, idx_gen, eopt)
%applyPTC: Apply production subsidy to certain generators

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Do not apply to DL
idx_gen = idx_gen & ~strcmp(mpc.genfuel, 'dl');

%% Apply credit
mpc.gen_aux{idx_gen, 'PTC'} = mpc.gen_aux{idx_gen, 'PTC'} + credit;
mpc = updateGenCost(mpc, idx_gen, eopt);

%% Display progress information
vfprintf(eopt.verbose, 'Production subsidy of $%.2f applied to gencost for %d generators\n', credit, sum(idx_gen));
