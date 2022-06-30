function offer = syncOffer(mpc, offer, eopt)
%syncOffer: Sync offer with mpc.gen
% 	For example, if dispatchable loads scaled.
%   DAC units use negative active reserve quantity, not positive, like
%   loads and other generators do.
%
% Inputs
%   MPC - Standard MATPOWER case struct
% 	OFFER - Standard offer matrix
% 	EOPT - E4ST Options Struct
%
% Examples
%   offer = syncOffer(mpc, offer, eopt);

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Sync offer
define_constants;
idx_dac = strcmp(mpc.genfuel, 'dac');
offer(~idx_dac, 2) = max(mpc.gen(~idx_dac, PMAX), -mpc.gen(~idx_dac, PMIN));
offer(idx_dac, 4) = abs(mpc.gen(idx_dac, PMIN));

vfprintf(eopt.verbose, 'Offer synced with PMAX\n');
