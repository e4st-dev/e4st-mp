function [mpc, esc] = addAreas(mpc, esc)
% addAreas Add custom areas to e4st bus and generator mappings
% - Custom areas defined in esc.custom_map
% - Updated mappings include esc.bus_map, mpc.gen_map

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Initialize
info = esc.custom_map;
info = filterInfo(info, esc);

if isfield(info, 'joint')
    info_list = info.joint{:, :};
else
    info_list = (1:height(info.value));
end

[bus_areas, gen_areas] = deal([]);

%% Add each new area to mappings
for i = unique(info_list)'
    idx_info = i == info_list;

    cur_info = filterStruct(info, idx_info);

    idx_bus = getInfoIdx(mpc, esc, cur_info, 'bus');
    idx_gen = getInfoIdx(mpc, esc, cur_info, 'gen');

    % Map idx_bus from mapping in mpc.bus to esc.bus_map
    % getInfoIdx uses mpc.bus to get index
    bus_idx = array2table([mpc.bus(:, 1), idx_bus], 'VariableNames', {'bus', 'idx'});
    bus_idx = join(esc.bus_map(:, 'bus'), bus_idx);
    idx_bus = bus_idx{:, 'idx'} == 1;

    new_area = unique(cur_info.value{:, 'new_area'});
    new_subarea = unique(cur_info.value{:, 'new_subarea'});
    new_default = unique(cur_info.value{:, 'new_default'});
    type = unique(cur_info.type{:, 'type'});

    % Bus mappings
    if any(strcmp(type, {'busgen', 'bus'}))
        if ~ismember(new_area, esc.bus_map.Properties.VariableNames) || ~any(strcmp(new_area, bus_areas))
            esc.bus_map{:, new_area} = new_default;
            bus_areas = [bus_areas; new_area];
        end
        esc.bus_map{idx_bus, new_area} = new_subarea;
    end
    
    % Generator mappings
    if any(strcmp(type, {'busgen', 'gen'}))
        if ~ismember(new_area, mpc.gen_map.Properties.VariableNames) || ~any(strcmp(new_area, gen_areas))
            mpc.gen_map{:, new_area} = new_default;
            gen_areas = [gen_areas; new_area];
        end
        mpc.gen_map{idx_gen, new_area} = new_subarea;
    end

end
