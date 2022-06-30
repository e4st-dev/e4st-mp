function res = ccs_retrofit_curves(var, avgcap, hr)
%CCS_RETROFIT_CURVES
%Function used to estimate cost of CCS retrofit and the heat rate penalty
%from the retrofit. These fuctions are derived from a regression in EPA
%Schedule 6 data. For more information, see the folder:
%
%Inputs:
%   Variable Name
%   Capacity: MW
%   Heat Rate (HR): mmbtu/mwh
%
%Outputs
%   Heat Rate Penalty: (%)
%   Capacity Penalty: (%)
%   Capital Cost: (2016$/KW)
%   FOM ($/kW-yr)
%   VOM ($/MWh)

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%%  Function code 
%truncate to avoid extrapolating EPA data. Power plants that are larger
%than the largest EPA class will simply have the same properties as the
%largest EPA class. 
idx = avgcap > 1000; avgcap(idx) = 1000;
idx2 = avgcap < 300; avgcap(idx2) = 300;

switch var
    case 'CAP_COST'
        res = 2496.4444 + -6.9022 * avgcap + 0.003544 * avgcap.^2 + 267.3333 * hr;
    case 'FOM'
        res =  48.503704 + -0.116685 * avgcap + 0.0000598148 * avgcap.^2 +  3.000000  * hr;
    case 'VOM'
        res = 1.4281 + -0.0060963 * avgcap + 0.0000031481 * avgcap.^2 + 0.4233333333 * hr;
    case 'CapPen' %capacity penalty
        res = 49.2333 + -0.112000 * avgcap + 0.0000533333 * avgcap.^2 + 2.4333333333 * hr;
    case 'HRPen'  %heat rate penalty of the retrofit
        res = 89.774 + -0.2513148 * avgcap + 0.00012907 * avgcap.^2 + 5.00000 * hr;
end
end

