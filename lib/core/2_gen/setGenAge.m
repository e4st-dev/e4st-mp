function mpc = setGenAge(mpc, year_delta, eopt)
% setGenAge: Update age of generators
%
% Inputs
%   MPC - Standard MATPOWER case struct
%   YEARDELTA - Number of years between current and previous simulation
% 	EOPT - E4ST Options structure

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Update generator age
mpc.gen_desc{:, 'AGE'} = mpc.gen_desc{:, 'AGE'} + year_delta;

%% Update generator vintage
mpc.newgen = mpc.newgen - 1;

%% Display progress information
vfprintf(eopt.verbose, 'Generators aged by %d years\n', year_delta);
