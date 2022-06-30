function [mpc, contab] = makePrlNew(mpc, contab, esc, eopt)
%makePrlNew: Create gencost steps for price-responsive load in base hour and contingency hours

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Extend mpc columns if necessary
% 10 steps of pairs; first pair is already included
if size(mpc.gencost, 2) ~= 24
    mpc.gencost = [mpc.gencost(:, 1:5), zeros(mpc.ng, 19)];
end

%% Make step gencost for base hour
idx_dl = strcmp(mpc.genfuel, 'dl');
mpc.gencost(idx_dl, :) = makePrlStep(mpc, esc, 1);

%% Make step gencost of contab for contingency hours
contab_prl = makePrlCont(mpc, esc);
contab = [contab; contab_prl];

%% Display progress information
vfprintf(eopt.verbose, 'Price responsive load set in year %d\n', esc.year);
