function mpc = matchLoad(mpc, esc, eopt)
%matchLoad: Scale loads to match provided values by location

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Prepare Load Match Data
define_constants;

info = esc.load_match;
info = filterInfo(info, esc);
info.value = info.value(:, ['Y', num2str(esc.year)]);

nh = size(esc.hrs_map, 1);
load_zone = unique(mpc.bus(:, BUS_AREA));
hourly_shape = mpc.hourly_shape;

%% Scale loads
% bus_base_load = total_load(mpc.bus, mpc.gen, (1: size(mpc.bus, 1))');
% bus_ann_load = bus_base_load .* hourly_shape * esc.hrs_map{:, 'hours'};

opt = struct('pq', 'P', 'scale', 'QUANTITY');
for i = 1:height(info.value)
    bus_base_load = total_load(mpc.bus, mpc.gen, (1:size(mpc.bus, 1))');
    bus_ann_load = bus_base_load .* hourly_shape * esc.hrs_map{:, 'hours'};

    idx_info = i == (1:height(info.value))';
    cur_info = filterStruct(info, idx_info);

    idx_bus = getInfoIdx(mpc, esc, cur_info, 'bus');

    if ~any(idx_bus)
        continue
    end

    tgt_ann_load = cur_info.value{:, :};

    zones = unique(mpc.bus(idx_bus, BUS_AREA));
    for zone = zones'
        idx_bus_zone = idx_bus & mpc.bus(:, BUS_AREA) == zone;
        zone_share = sum(bus_ann_load(idx_bus_zone)) / sum(bus_ann_load(idx_bus));
        zone_tgt = tgt_ann_load * zone_share;

        if zone_share == 0
            continue
        end
        cur_base_load = sum(bus_base_load(idx_bus_zone));
        cur_ann_load = sum(bus_ann_load(idx_bus_zone));
        ratio = cur_base_load / cur_ann_load;
        tgt_base_load = zone_tgt * ratio;

        [mpc.bus, mpc.gen] = scale_load(tgt_base_load, mpc.bus, mpc.gen, idx_bus_zone, opt);
    end
    vfprintf(eopt.verbose, 'Load scaled in %d buses to match %.2f annual load\n', sum(idx_bus), tgt_ann_load);
end

%% Calculate total annual loads by bus
hourly_load = zeros(size(mpc.bus, 1), nh);
ann_load = zeros(size(mpc.bus, 1), 1);
base_load = total_load(mpc.bus, mpc.gen, (1:size(mpc.bus, 1))');
for i_area = 1:length(load_zone)
    idx_bus = mpc.bus(:, BUS_AREA) == i_area;
    hourly_load(idx_bus, :) = base_load(idx_bus) .* hourly_shape(idx_bus, :);
    ann_load(idx_bus) = hourly_load(idx_bus, :) * esc.hrs_map{:, 'hours'};
end

mpc.hourly_load = hourly_load;
mpc.ann_load = ann_load;
