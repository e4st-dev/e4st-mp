function mpc = updateMPCDim(mpc)
% updateMPCDim Update dimensions of MPC total_output and caplim constraints

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Update bus and generator size
mpc.nb = size(mpc.bus, 1);
mpc.ng = size(mpc.gen, 1);

%% Update constraint dimensions
% total output constraint
if isfield(mpc, 'total_output') && isfield(mpc.total_output, 'map') && ...
   ~isempty(mpc.total_output.map) && size(mpc.total_output.map, 2) < mpc.ng
    map = mpc.total_output.map;
    mpc.total_output.map = [map  ...
                            zeros(size(map, 1), mpc.ng - size(map, 2))];

    mpc.total_output.coeff = [mpc.total_output.coeff;  ...
                            zeros(mpc.ng - size(map, 2), size(map, 1))];
end

% capacity limit constraint
if isfield(mpc, 'caplim') && isfield(mpc.caplim, 'map') && ...
   ~isempty(mpc.caplim.map) && size(mpc.caplim.map, 2) < mpc.ng
    map = mpc.caplim.map;
    mpc.caplim.map = [map  ...
                            zeros(size(map, 1), mpc.ng-size(map, 2))];
end
