function [mpc, contab] = create_contab(mpc, esc, eopt)
% create_contab creates the contingency table and modifies the MPC
% file for use in the MATPOWER superopf. 
%As input, it takes the standard e4st case file (esc) and
%matpower case file (mpc). It also reads in the e4st options structure
%(eopt), in case it is useful for future function modifications.
%The function works by reading hour probabilities from the esc file and
%hourly loads at each bus from the mpc file (mpc.hourly_load). It then
%converts the hourly_load into a contingency table, known internally as
%'contab'. This contingency table is essentially a detailed list of
%changes that must be applied to the model form one representative hour to
%the next representative hour. For more details on the contingency table,
%type help apply_changes, or help idx_ct in a MATPOWER activated terminal.
%In the mpc file, we change the mpc.bus(:, BUS_AREA) and the base generator
%loads.

% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Christoph Funke (Resources for the Future)

% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).


%% Begin function

vfprintf(eopt.verbose, "Creating contingency table ... ")

%% Setup
define_constants;

% Type 'help apply_changes' and 'help idx_ct' for details on change table
nh = size(esc.hrs_map, 1); %number of hours
nb = height(mpc.bus); %number of buses

%hour probabilities
hr_prob = esc.hrs_map{:, 'probability'};
tol = 10^-14;
assert(abs(sum(hr_prob) - 1) < tol, 'Hour probabilities do not sum to 1')
mpc.hr_prob = hr_prob;

%innitialize change table
contab = zeros(nb*(nh-1), 7);
mpc.bus(:, BUS_AREA) = [1:nb]; %give each bus it's own area

%% Scale base hour to adjust for added loads

opt = struct('pq', 'P');
old_base_load = total_load(mpc.bus, mpc.gen, (1:size(mpc.bus, 1))');
new_base_load = mpc.hourly_load(:,1);
tol= 10^-10;
idx_update = abs(new_base_load - old_base_load) > tol; %update only if load changed

for i = find(idx_update)'
    if ~any(old_base_load(i)) %skip if no load at bus or no buses to be updated
        continue
    end
    idx_bus = i == (1:nb)';
    load_factor = new_base_load(i)./old_base_load(i);
    [mpc.bus, mpc.gen] = scale_load(load_factor, mpc.bus, mpc.gen, idx_bus, opt);
    vfprintf(eopt.verbose, 'Base loads scaled by %.4f at %d bus to acommodate added load\n', load_factor, sum(idx_bus));
end

%check if successful
new_load = total_load(mpc.bus, mpc.gen, (1:size(mpc.bus, 1))');
assert(all(abs(new_base_load - new_load) < tol, 'all'), "Baseload set incorrectly in create_contab") ;

%% Create contingency table 
%for a desctiption of the contab format, type 'help apply_changes' and 'help
%idx_ct' in a MATPOWER activated terminal window
i_contab = 1;
for i_bus = 1:nb
    bus_area = i_bus; % since bus area was defined above as being equal to the row index in mpc.bus
        
    for i_hr = 2:nh
        hr_label = i_hr - 1;
        prob = hr_prob(i_hr); 
        tableID = 8;  %area-wide load changes (CT_TAREALOAD in Table 9-2)
        col = 4;  %modify all loads, real only (CT_LOAD_ALL_P  in idx_ct)
        type = 1; %replaces old value by value in CT_NEWVAL column (CT_REP in idx_ct)
        value = mpc.hourly_load(i_bus, i_hr); %need to break this into hourly values
        contab(i_contab,:) = ...
            [hr_label, prob, tableID, bus_area, col, type, value];
        i_contab = i_contab + 1;
    end
end


%% Print Results

vfprintf(eopt.verbose, " done\n");
end

