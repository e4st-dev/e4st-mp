function mpc = applyIRMReq(mpc, pol_val, idx_gen, pol_name, eopt)
%applyIRMReq: apply Installed Reserve Margin policies

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Set Installed Reserve Requirement (Capacity Constraint)
%define_constants;
PMAX = 9;
[map, capMax, capMin] = deal([]);

%%
if sum(mpc.gen(:, PMAX).*idx_gen) >= pol_val
    map = [map; idx_gen'];
    capMax = [capMax; Inf];
    capMin = [capMin; pol_val];
    vfprintf(eopt.verbose, 'IRM capacity requirement set to %.2f MW for %d generators\n', pol_val, sum(idx_gen > 0));
    mpc = addCapLim(mpc, map, capMin, capMax, pol_name);
else
    fprintf('Warning: Not enough capacity to meet IRM requirement\n')
    fprintf('\t Requirement: %.2f MW\n', pol_val)
    fprintf('\t Eligible Capacity: %.2f\n', sum(mpc.gen(:, 2) .* idx_gen))
    %error('Not enough capacity to meet IRM requirement\n')
end
