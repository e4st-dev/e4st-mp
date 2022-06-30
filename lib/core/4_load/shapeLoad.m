function mpc = shapeLoad(mpc, esc, eopt)
%shapeLoad: Set load shape by bus and write it to mpc.hourly_load
%
% Inputs
%   MPC - Standard MATPOWER case struct
%
%   ESC - Standard E4ST case struct
%
%   eopt - Standard e4st options file
%
% Examples
%		mpc = shapeLoad(mpc, esc, eopt)

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano, Christoph 
% Funke, and Steven Witkin (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Prepare Load Shape Data
info = esc.load_shape;
info = filterInfo(info, esc);

nh = size(esc.hrs_map, 1);

%% Update Bus Area
define_constants;

%% Calculate total annual loads by bus
hourly_shape = zeros(size(mpc.bus, 1), nh);

%identify buses in each region of load_shape and assign each region a
%separate BUS_AREA
for i = 1:height(info.value)
    idx_info = i == (1:height(info.value))';
    cur_info = filterStruct(info, idx_info);
    idx_bus = getInfoIdx(mpc, esc, cur_info, 'bus');
    
    area = cur_info.area{:, 'area'};
    subarea = cur_info.area{:, 'subarea'};  
    values = cur_info.value{:,:};
    
    
    if ~any(idx_bus)
        continue
    end
    hourly_shape(idx_bus, :) = repmat(values, sum(idx_bus), 1);
    
    vfprintf(eopt.verbose, "Load shape set for %d buses in %s %s\n", sum(idx_bus), area{1}, subarea{1});
end

mpc.hourly_shape = hourly_shape;

%% Display progress information
vfprintf(eopt.verbose, 'Scale loads in %d area(s) in %d contingency hours\n', height(info.value), nh-1);
