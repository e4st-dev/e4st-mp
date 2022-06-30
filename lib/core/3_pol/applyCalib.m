function mpc = applyCalib(mpc, esc, calib_val, calib_type, idx_gen, pol_name, eopt)
%applyCalib: Apply price or quantity calibration policy

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Apply Calibration
%define_constants;
if strcmp(calib_type, 'prc') || strcmp(calib_type, 'price')
    
    % Standard variable cost calibration
    mpc.gen_aux{idx_gen, 'PTC'} = mpc.gen_aux{idx_gen, 'PTC'} - calib_val;
    mpc = updateGenCost(mpc, idx_gen, eopt);
    vfprintf(eopt.verbose, 'Calibration cost adder of %.3f applied to %d generators\n', calib_val, sum(idx_gen));
    
elseif strcmp(calib_type, 'qty')
    cap = calib_val / 8760; % hourly output requirement
    coeff = ones(length(idx_gen), 1) .* idx_gen;

    max_output = sum(mpc.gen(:, 9).*coeff.*(mpc.availability_factor(:, :) * esc.hrs_map{:, 'hours'})) / 8760; %PMAX = 9
    if max_output < cap % Check if feasible
        vfprintf(eopt.verbose, 'Calibration constraint of %.3f is infeasible, switching to max output of %.3f\n', calib_val, max_output);
        cap = max_output;
    end

    map = idx_gen';
    min = cap * .99; % -Inf;
    max = cap * 1.01;
    type = 1;
    name = pol_name;

    mpc = addTOC(mpc, map, min, max, coeff, type, name);
    %mpc = addTOC(mpc, idx_gen', cap*-Inf, cap, coeff, 1, {pol_name});
    vfprintf(eopt.verbose, 'Calibration constraint of %.3f applied to %d generators\n', calib_val, sum(idx_gen));
end
