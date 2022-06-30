function offer = updateOfferPrc(mpc, offer, idx_gen, eopt)
% updateOfferPrc Update offer price ($/MW/hr) for specified generators
% - offer price = sum(FOM, CAP_COST, TRANS_COST, ROUTINE_CAPEX) - ITC
% - idx_gen specifies indices for generators to update

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Check generator index
if isempty(idx_gen)
    idx_gen = ones(mpc.ng, 1) == 1;
end

%% Do not update dispatchable load, and treat DAC differently
idx_gen = idx_gen & ~strcmp(mpc.genfuel,'dl');
idx_dac = idx_gen & strcmp(mpc.genfuel, 'dac');
idx_gen = idx_gen & ~idx_dac;

%% Update offer price: offer price = sum(FOM, CAP_COST, TRANS_COST, ROUTINE_CAPEX) - ITC
offer(idx_gen, 1) = sum(mpc.gen_aux{idx_gen, {'FOM', 'CAP_COST', 'TRANS_COST', 'ROUTINE_CAPEX'}}, 2) ...
                    - mpc.gen_aux{idx_gen, 'ITC'};
                
offer(idx_dac, 3) = sum(mpc.gen_aux{idx_dac, {'FOM', 'CAP_COST', 'TRANS_COST', 'ROUTINE_CAPEX'}}, 2) ...
                    - mpc.gen_aux{idx_dac, 'ITC'};

%% Display progress information
vfprintf(eopt.verbose, '(Offer price updated for %d generators)\n', sum(idx_gen) + sum(idx_dac));
