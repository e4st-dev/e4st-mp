function [mpc, esc] = selIsland(mpc, esc, island, eopt)
% selIsland Select islands from MPC to eliminate isolated buses

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Extract islands
define_constants;
if length(island) == 1 && island == 0
    island = 'all';
end

% Temporary index of mpc gen tables (extract islands does not accept tables)
mpc.gen_idx = (1:mpc.ng)';
% Temporary index of mpc gen tables (extract islands does not accept tables)
mpc.bus_idx = (1:mpc.nb)';

nb = mpc.nb;

% Extract islands
precap = sum(mpc.gen(:, PMAX));
custom.bus{1} = {'bus_island', 'bus_name', 'bus_idx'};
custom.gen{1} = {'newgen', 'gen_idx'};
mpc = extract_islands(mpc, island, custom);
postcap = sum(mpc.gen(:, PMAX));

% Check if generaiton capacity removed
if precap ~= postcap
    error('Generation capacity located on island and removed\n\tPre-Cap: %.2f, Post-Cap: %.2f, Removed Cap: %.2f\n', precap, postcap, postcap-precap);
end

% Filter gen tables
gen_tables = {'gen_aux', 'gen_map', 'gen_desc'};
for field = gen_tables
    mpc.(char(field)) = mpc.(char(field))(mpc.gen_idx, :);
end
mpc = rmfield(mpc, 'gen_idx');

% Filter bus tables
if isfield(esc, 'bus_map')
    esc.bus_map = esc.bus_map(mpc.bus_idx, :);
end
mpc = rmfield(mpc, 'bus_idx');

mpc.nb = size(mpc.bus, 1);
mpc.ng = size(mpc.gen, 1);
nRm = nb - mpc.nb;

%% Display progress information
if ischar(island)
    vfprintf(eopt.verbose, '%s islands (%d buses) extracted and removed from MPC\n', island, nRm);
else
    vfprintf(eopt.verbose, '%d islands (%d buses) extracted and removed from MPC\n', island, nRm);
end
