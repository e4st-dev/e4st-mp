function [mpc, offer] = applyITC(mpc, offer, credit, idx_gen, eopt)
%applyITC: Apply investment subsidy to certain generators

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Do not apply to DL
idx_gen = idx_gen & ~strcmp(mpc.genfuel, 'dl');

%% Apply credit
mpc.gen_aux{idx_gen, 'ITC'} = mpc.gen_aux{idx_gen, 'ITC'} + mpc.gen_aux{idx_gen, 'CAP_COST'} * (credit);
offer = updateOfferPrc(mpc, offer, idx_gen, eopt);

%% Display progress information
vfprintf(eopt.verbose, 'Investment tax credit of %.2f percent applied to %d generators\n', credit*100, sum(idx_gen));
