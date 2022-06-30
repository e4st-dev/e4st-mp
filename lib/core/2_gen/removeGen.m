function mpc = removeGen(mpc, idx_gen)
%removeGen: Remove a set of generators from mpc and offer; Update constraint dimensions
%
% Inputs
%   MPC - Standard MATPOWER case struct
%   idx_gen - Index of generators to remove

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%global verbose % Display Progress Information
verbose = 1;

%% Update Main MPC Fields
mpc_fields = {'gen', 'gencost', 'genfuel', 'gentype', 'newgen', ...
    'gen_aux', 'gen_desc', 'gen_map', 'availability_factor', 'short_term_storage'};

for field = mpc_fields
    if isfield(mpc, field)
        mpc.(char(field))(idx_gen, :) = [];
    end
end

mpc.ng = size(mpc.gen, 1);

%%

% Update Emission Controls
if isfield(mpc, 'emis_ctrl')
    if ~isempty(mpc.emis_ctrl)
        for field = fieldnames(mpc.emis_ctrl)
            mpc.emis_ctrl.(char(field))(idx_gen, :) = [];
        end
    end
end

% Update Total_Output Constraints
if isfield(mpc, 'total_output')
    if isfield(mpc.total_output, 'map')
        if ~isempty(mpc.total_output.map)
            mpc.total_output.map(:, idx_gen) = [];
            mpc.total_output.coeff(idx_gen, :) = [];
        end
    end
end

% Update Capacity Constraints
if isfield(mpc, 'caplim')
    if isfield(mpc.caplim, 'map')
        if ~isempty(mpc.caplim.map)
            mpc.caplim.map(:, idx_gen) = [];
        end
    end
end

%% Display Progress Information
vfprintf(verbose, 'A set of %d generators has been removed\n', sum(idx_gen));
