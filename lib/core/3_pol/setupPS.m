function [idx_gen, idx_dl, pol_val] = setupPS(mpc, esc, idx_gen, gen_wgt, cur_info, pol_val)
%setupPS prepares RPS/CES policy specification

%Inputs
%idx_gen is index of generators relevant to policies.
%gen_wgt is maximum capacity of generators that can count toward policies.
%cur_info is the data rows in esc.policy that is relevant for the policy.
%pol_val = the value of the policy.

%returns index of relevant generators (idx_gen)
%returns list with policy_value at generators of type 'dl'(idx_dl)
%changes policy_value to 1

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

% if any indices in idx_gen refer to demand side (dl, DAC, storage)
if any(idx_gen & ismember(mpc.genfuel, {'dl'}))
    % Remove things purely on demand side from idx_gen
    idx_gen = idx_gen & ~(ismember(mpc.genfuel, {'dl', 'storage'})); %| strcmp(mpc.genfuel, 'dac')); 
    % idx_dl is index of demand side
    idx_dl = gen_wgt .* ismember(mpc.genfuel, {'dl', 'dac', 'storage'});
    pol_val = 1; %sum(gen_wgt(idx_dl>0) .* mpc.gen(idx_dl>0, 2)) / sum(mpc.gen(idx_dl>0, 2));
else
    %index of generators in area
    idx_area_gen = getInfoIdx(mpc, esc, rmfield(cur_info, {'genfuel', 'gentype'}), 'gen');
    %note mpc.ng = number of generators :)
    %idx_dl returns pol_val if generator type 'dl' is in the area
    idx_dl = pol_val * ones(mpc.ng, 1) .* idx_area_gen .* ismember(mpc.genfuel, {'dl', 'dac', 'storage'}); %(strcmp(mpc.gentype, 'dl') | strcmp(mpc.genfuel, 'dac'));
    idx_stor = strcmp(mpc.genfuel, 'storage');
    idx_dl(idx_stor) = idx_dl(idx_stor) * -0.15/0.85;
    idx_gen = idx_gen | idx_area_gen .* strcmp(mpc.genfuel, 'dac');
    pol_val = 1;
end

%convert policy value to actual value
if ~(pol_val >= 0 && pol_val <= 1) % target value, convert to fraction
    idx_area_bus = getInfoIdx(mpc, esc, rmfield(cur_info, {'genfuel', 'gentype'}), 'bus');
    pol_val = pol_val / sum(mpc.ann_load(idx_area_bus));
end
end