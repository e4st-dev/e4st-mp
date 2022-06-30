function mpc = addCapLim(mpc, map, min, max, cap_name)
% addCapLim Add capacity limit (minimum and/or maximum)

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Paul Picciano (Resources for the Future) and Biao Mao (Rensselaer Polytechnic Institute)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Check if caplim field is initialized
if isempty(mpc.caplim) && ~isempty(map)
    mpc.caplim.map = [];
    mpc.caplim.max = [];
    mpc.caplim.min = [];
    mpc.caplim.name = [];
end

%% Check if caplim map is correct size
% If new generators have been added recently, old map is too small
[a,b] = size(mpc.caplim.map);
if a > 0
    if b < length(map)
        mpc.caplim.map = [mpc.caplim.map, zeros(a, length(map)-b)];
    elseif b == length(map)
        %continue
    else    
        error("There is some issue. Generators deleted since last used caplim?")
    end    
end

%% Check name format
if ischar(cap_name)
    cap_name = {cap_name};
end

%% Update caplim constraints
if ~isempty(map)
    mpc.caplim.map = [mpc.caplim.map; map];
    mpc.caplim.min = [mpc.caplim.min; min];
    mpc.caplim.max = [mpc.caplim.max; max];
    mpc.caplim.name = [mpc.caplim.name; cap_name];
end
