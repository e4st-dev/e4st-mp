function loc_map = getMap(mpc, esc, map_type, map_area)
% getMap Return map for buses, generators, or branches
% - map_type specifies mapping (bus, gen)
% - map_area specifies area within mapping (if map_area = map_type or non
%   provided, full mapping returned)

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Get Mapping
% Generator map
if strcmp(map_type, 'gen')
    loc_map = table(mpc.gen(:, 1), 'VariableNames', {'bus'});
    if strcmp(map_area, 'bus')
    elseif strcmp(map_area, 'all_bus')
        loc_map{:, 'all_bus'} = {'in'};
    elseif isempty(map_area) && isfield(mpc, 'gen_map')
        loc_map = join(loc_map, esc.bus_map);
        for col = mpc.gen_map.Properties.VariableNames
            if any(strcmp(loc_map.Properties.VariableNames, col))
                loc_map(:, col) = [];
            end
            loc_map = [loc_map mpc.gen_map(:, col)];
        end
    elseif isempty(map_area)
        loc_map = join(loc_map, esc.bus_map);
    elseif isfield(mpc, 'gen_map') && ismember(map_area, mpc.gen_map.Properties.VariableNames)
        loc_map = mpc.gen_map(:, map_area);
    elseif ismember(map_area, esc.bus_map.Properties.VariableNames)
        loc_map = join(loc_map, esc.bus_map(:, {'bus', char(map_area)}));
    else
        error('Map area "%s" not recognized\n', char(map_area))
    end
% Bus map
elseif strcmp(map_type, 'bus')
    loc_map = table(mpc.bus(:, 1), 'VariableNames', {'bus'});
    if strcmp(map_area, 'bus')
    elseif strcmp(map_area, 'all_bus')
        loc_map{:, 'all_bus'} = {'in'};
    elseif isempty(map_area)
        loc_map = join(loc_map, esc.bus_map);
    elseif ismember(map_area, esc.bus_map.Properties.VariableNames)
        loc_map = join(loc_map, esc.bus_map(:, {'bus', char(map_area)}));
    else
        error('Map area "%s" not recognized\n', char(map_area))
    end
else
    error('Map type "%s" not recognized\n', char(map_type))
end
