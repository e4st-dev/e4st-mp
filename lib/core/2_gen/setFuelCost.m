function mpc = setFuelCost(mpc, esc, eopt)
%setFuelCost: Update generators fuel costs by fuel type

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Prepare Fuel Data
info = esc.fuel_cost;
info = filterInfo(info, esc);

idx_new = mpc.newgen == 1;

% some DAC pay for NG for heating. to achieve this without changing rows in
% esc, duplicate NG rows for DAC 
if strcmp(eopt.dac, 'T')
    idx_ng = strcmp(info.genfuel{:,:}, {'ng'});
    fields = fieldnames(info);
    for i = 1:length(fields)
        tmp = info.(char(fields(i)));
        if strcmp(fields(i), 'genfuel')
            tmp{idx_ng, 'genfuel'} = {'dac'};
        end
        info.(char(fields(i))) = [info.(char(fields(i))) ; tmp(idx_ng, :)];
    end
end

%% Set Fuel Data
has_fuel = zeros(mpc.ng, 1); %variable to track whether fuel cost has been updated for a given generator

for i = 1:height(info.value)
    idx_info = i == (1:height(info.value))';
    cur_info = filterStruct(info, idx_info);
    idx_gen = getInfoIdx(mpc, esc, cur_info, 'gen');

    if ~any(idx_gen)
        continue
    end
   
    % Pre-Existing
    mpc.gen_aux{idx_gen & ~idx_new, 'FUEL_COST'} = ...
        mpc.gen_aux{idx_gen & ~idx_new, 'FUEL_COST'} + ...
        cur_info.value{:, ['Y', num2str(esc.year)]} - cur_info.value{:, ['Y', num2str(esc.year - esc.year_delta)]};

    % New Generators
    mpc.gen_aux{idx_gen & idx_new, 'FUEL_COST'} = cur_info.value{:, ['Y', num2str(esc.year)]};

    % Update Gencost
    mpc = updateGenCost(mpc, idx_gen, eopt);

    vfprintf(eopt.verbose, 'Fuel cost and gencost updated for %d %s generators\n', sum(idx_gen), string(cur_info.genfuel{:,:}));
    has_fuel = has_fuel | idx_gen; %mark generator as being updated
end

%% check that all generators that require a fuel cost have been updated
% Assume that all online generators with a non-zero heat rate must have a fuel
% cost
tol = 10^-5;
idx = ~has_fuel & mpc.gen_aux{:, 'HR'} > tol & mpc.gen(:,8); 
if sum(idx) ~= 0
    genfuel = mpc.genfuel{find(idx,1)};
    error(['No fuel cost assigned to generator with fuel type: ' genfuel])
end