function mpc = addTOC(mpc, map, min, max, coeff, type, toc_name)
% addTOC Add total output constraint (minimum and/or maximum)
% - formulation: sum(generator output * coeff) [< >] K
% - map: mapping to apply TOC (optional third dimension is representative hour)
% - min: minimum limit
% - max: maximum limit
% - coeff: coefficient to multiply generator output
% - type: column
% - name: name of applied toc for e4st results processing

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Paul Picciano (Resources for the Future) and Biao Mao (Rensselaer Polytechnic Institute)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Check if total_output filed is initialized
if isempty(mpc.total_output)
    mpc.total_output.map = [];
    mpc.total_output.max = [];
    mpc.total_output.min = [];
    mpc.total_output.coeff = [];
    mpc.total_output.type = [];
    mpc.total_output.name = [];
end

%% Check if map and coeff are correct size
% If new generators have been added recently, old map is too small
[a,b] = size(mpc.total_output.map);
if a > 0
    if b < length(map)
        mpc.total_output.map = [mpc.total_output.map, zeros(a, length(map)-b)];
        mpc.total_output.coeff = [mpc.total_output.coeff; zeros(length(map)-b, a)];
    else
        % there is some issue. generators deleted since last used TOC?
    end    
end

%% Check name format
if ischar(toc_name)
    toc_name = {toc_name};
end

%% Hourly map
% Warning: Substantially increases RAM usage
% if size(map, 3) == 1
%     % Set map in all hours
%     map = repmat(map, 1,1,mpc.nh);
% elseif size(map, 3) ~= mpc.nh
%     % Incorrect hours, clear map to avoid applying constraint
%     map = [];
% end

%% Update total_output constraints
if ~isempty(map)
    mpc.total_output.map = [mpc.total_output.map; map];
    mpc.total_output.min = [mpc.total_output.min; min];
    mpc.total_output.max = [mpc.total_output.max; max];
    mpc.total_output.coeff = [mpc.total_output.coeff coeff];
	mpc.total_output.type = [mpc.total_output.type; type * size(mpc.total_output.coeff, 2)];
    mpc.total_output.name = [mpc.total_output.name; toc_name];
end
