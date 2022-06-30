function offer = setGenOffer(mpc, eopt)
%setGenOffer: Create offer table and set offers for all pre-existing generators
%
%       offer
%           .PositiveActiveReservePrice
%           .PositiveActiveReserveQuantity
%           .NegativeActiveReservePrice
%           .NegativeActiveReserveQuantity
%           .PositiveActiveDeltaPrice
%           .NegativeActiveDeltaPrice
%           .PositiveActiveReservePrice2        (optional quadratic term)
%           .NegativeActiveReservePrice2        (optional quadratic term)
%           .PositiveActiveDeltaPrice2          (optional quadratic term)
%           .NegativeActiveDeltaPrice2          (optional quadratic term)
%           .ActiveContractMin                  (optional)
%           .ActiveContractMax                  (optional)
%           .PminFactor                         (optional)
%
% Inputs
%   MPC - Standard MATPOWER case struct
% 	EOPT - E4ST Options Struct

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Intialize the offer table
offer = zeros(mpc.ng, 13);
offer(:, 4) = Inf;
%offer(:, 11) = -Inf; % contract pc_min
%offer(:, 12) = Inf; % contract pc_max

%% Set offer for all fuel types
% Direct air capture is different, in that its offer is [0, 0, offerPrc, buildableCap]
% rather than other gens, which have [offerPrc, buildableCap, 0, Inf]
define_constants;
% Set the PositiveActiveReservePrice % Fixed Costs, including tax, insurance and routine capex
offer = updateOfferPrc(mpc, offer, [], eopt);
offer(strcmp(mpc.genfuel, 'dl'), 1) = 0;

% Set the PositiveActiveReserveQuantity for non-dac, and NegativeARP for DAC
idx_dac = strcmp(mpc.genfuel, 'dac');
offer(~idx_dac, 2) = max(mpc.gen(~idx_dac, PMAX), -mpc.gen(~idx_dac, PMIN));
offer(idx_dac, 4) = abs(mpc.gen(idx_dac, PMIN));

%% Display progress information
vfprintf(eopt.verbose, 'Offer set for all pre-existing generators \n');
