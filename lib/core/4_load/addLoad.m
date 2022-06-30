function [mpc, contab] = addLoad(mpc, esc, eopt)
%addLoad add load at specific bus/hours

% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Set UP
define_constants;

info = esc.load_add;
info = filterInfo(info, esc); %get only relevant rows (filter for sim year and status)

nh = size(esc.hrs_map, 1);
nb = size(mpc.bus,1);
hr_prob = esc.hrs_map{:, 'probability'};
tol = 10^-14;
assert(abs(sum(hr_prob) - 1) < tol, 'Hour probabilities do not sum to 1')


%% Create hourly_add table
%base load at each bus;
bus_base_load = total_load(mpc.bus, mpc.gen, (1:size(mpc.bus, 1))');


mpc.hourly_add = zeros(nb, nh); %initialize matrix to store added values;
for i_area = 1:height(info.area)
%    idx_info = i_area == (1:height(info.value))';
    cur_info = filterStruct(info, i_area);
    idx_bus = getInfoIdx(mpc, esc,  cur_info, 'bus');
    area = cur_info.area.area{1};
    % if subarea is a double, the cell indexing won't work
    try
        subarea = cur_info.area.subarea{1};
    catch
        subarea = num2str(cur_info.area.subarea);
    end
    nbus = sum(idx_bus);
    value = cur_info.value{:,:};
    
    base_loads = bus_base_load(idx_bus);
    assert(sum(base_loads) ~= 0, 'No current load in load_add subarea %s', subarea)
    bus_weights = base_loads / sum(base_loads);
    add_values = bus_weights * value;
    mpc.hourly_add(idx_bus, :) = mpc.hourly_add(idx_bus, :) + add_values;
      

    vfprintf(eopt.verbose, 'Load added in %s %s with %d buses\n', area, subarea, nbus);

end

%% Update mpc.hourly_shape, mpc.hourly_load, and mpc.ann_load
%these are importantly used in the results processing .... 

%hourly load
mpc.hourly_load = mpc.hourly_load + mpc.hourly_add; %add to hourly load;
assert(all(mpc.hourly_load >= 0, 'all'), 'Too much load subtracted in load_add. Negative or NA load values encountered.')

%hourly shape
mpc.old_hourly_shape = mpc.hourly_shape;
idx_update = any(mpc.hourly_load, 2); %ignore buses with no load
mpc.hourly_shape(idx_update,:) = mpc.hourly_load(idx_update,:) ./ mpc.hourly_load(idx_update, 1);
assert(~any(isnan(mpc.hourly_shape), 'all'), "Non-existent or missing value in mpc.hourly_shape")
assert(~nb ~= nh, "Hourly shape may have been incorrectly calculated. Check the calculation in addLoad")

%annual load
mpc.ann_load = mpc.hourly_load * esc.hrs_map{:, 'hours'};


