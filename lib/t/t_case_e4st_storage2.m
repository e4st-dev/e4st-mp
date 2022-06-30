function [mpc, roffer, contab] = t_case_e4st_storage2
%T_CASE2_E4ST_STORAGE   Two-bus test case for t_e4st_storage()
% Two-bus test case for testing e4st_solve() with short term storage.
% Contains two buses, each with a generator. The gererator at bus 1 has a 
% maximum power output of 10 MW and variable cost of $50/MW. The generator
% at bus 2 has a maximum power output of 20 MW and variable cost of $20/MW.
% There are two days. One day has 4 hours. For day 1, the first and third
% hour have a load of 15 MW, while the second and fourth hour have a load
% of 25 MW. For day 2, the first hour has a load of 25 MW and the second
% has a load of 15 MW. This test is designed to test whether days can have
% variable number of hours.


%% define named indices into data matrices
define_constants
c = idx_dcline;

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	0	0	0	1	1	0	135	1	1.05	0.95;
	2	1	0	0	0	0	2	1	0	135	1	1.1	0.95;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	10	0	5	-5	1	100	1	11	0	0	0	0	0	0	0	0	Inf	0	0	0;
	2	20	0	10	-10	1	100	1	21	0	0	0	0	0	0	0	0	Inf	0	0	0;
	% battery
 	2	10	0	0	0	1	100	1	11	-11/0.85	0	0	0	0	0	0	0	Inf	0	0	0;
	% dispatchable loads
	1	-10	0	0	0	1	100	1	0	-10	0	0	0	0	0	0	Inf	Inf	Inf	Inf	0;
	2	-5	0	0	0	1	100	1	0	-5	0	0	0	0	0	0	Inf	Inf	Inf	Inf	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.01	0.01	0.014	130	130	130	1	0	1	-360	360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	0	0	2	50	0;
	2	0	0	2	20	0;
	% battery
 	2	0	0	2	0	0;
	% dispatchable loads
	2	0	0	2	5000	0;
	2	0	0	2	5000	0;
];

%%-----  Short Term Storage Data  -----%%
%	gen_index	efficiency	max_energy(MWh/MW Capacity)
mpc.short_term_storage = [
	3	0.85	4;
];

%%-----  Day Data  -----%%
% each day is represented by a col vector
% The vector contains the labels of the contingencies that belong to that day,
% in the order in which they occur. Label 0 indicates base case.
mpc.days = {
	[0; 1; 2; 3], 24;
	[4; 5], 4;
};



%% reserve and delta offers
%     +      +        -      -
%  active  active  active  active
% reserve  reserve reserve reserve
%  price    qty    price    qty
roffer = [
    % generators
    10    10      1e-7    Inf;
    5     20      1e-7    Inf;
    % battery
    5      10      0       Inf;
%     % dispatchable loads
    1e-7    10      1e-7    Inf;
    1e-7    5       1e-7    Inf;
];

%% changes table
% label probty  field           row column          chgtype newvalue
contab = [
    1   0.125    CT_TAREALOAD    1   CT_LOAD_ALL_P   CT_REL  2;
    1   0.125    CT_TAREALOAD    2   CT_LOAD_ALL_P   CT_REL  1;
    2   0.125    CT_TAREALOAD    1   CT_LOAD_ALL_P   CT_REL  1.2;
    2   0.125    CT_TAREALOAD    2   CT_LOAD_ALL_P   CT_REL  1;
    3   0.125    CT_TAREALOAD    1   CT_LOAD_ALL_P   CT_REL  5;
    3   0.125    CT_TAREALOAD    2   CT_LOAD_ALL_P   CT_REL  1;
    4   0.25    CT_TAREALOAD    1   CT_LOAD_ALL_P   CT_REL  1.6;
    4   0.25    CT_TAREALOAD    2   CT_LOAD_ALL_P   CT_REL  1.2;
    5   0.25    CT_TAREALOAD    1   CT_LOAD_ALL_P   CT_REL  0;%same as original hour
    5   0.25    CT_TAREALOAD    2   CT_LOAD_ALL_P   CT_REL  0.8;
];

