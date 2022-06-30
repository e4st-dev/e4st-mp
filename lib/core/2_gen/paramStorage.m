function [mpc, offer] = paramStorage(mpc, esc, offer, eopt)
% paramStorage set up diurnal storage units

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Apply Storage Charging Parameters
%Set relevant column indices
define_constants;
GID = 1; %Generator Indices
EFF = 2; %Storage Efficiency
SCAP = 3; %Storage Energy Capcity
ng = size(mpc.gen,1); %number of generators

%Actual code
idx_storage = strcmp(mpc.genfuel(:,:), 'storage') & (mpc.gen(:,GEN_STATUS) == 1);
mpc.gen(idx_storage, PMIN) = -mpc.gen(idx_storage, PMAX)./...
    mpc.short_term_storage(idx_storage, EFF); %set pMin = -PMAX/eff
mpc.gen_aux{idx_storage, {'FOM', 'CAP_COST'}} = mpc.gen_aux{idx_storage, {'FOM', 'CAP_COST'}}.*...
    mpc.short_term_storage(idx_storage, EFF); %adjsut capital cost so that it is based on charging capacity not discharging capacity
offer = updateOfferPrc(mpc, offer, idx_storage, eopt);

mpc.short_term_storage(:, GID) = 1:ng;%reset generator index properly
mpc.short_term_storage = mpc.short_term_storage(idx_storage,:); %keep just storage units


%check if indexes are correct
check_storage = all(strcmp(mpc.genfuel(mpc.short_term_storage(:,GID)',:), 'storage'));
assert(check_storage, 'Incorrect generator indices in mpc.short_term_storage')

%% Group RprHrs into Days
if size(mpc.short_term_storage,1) ~=0
    hr_day = esc.hrs_map{:, 'e4st_day'};
    nd = size(unique(hr_day),1);
    mpc.days = cell(nd,2);
    for day = 1:nd
        %minus one makes the day numbers match the contingencies in contab
        idx = find(hr_day == day);
        mpc.days{day, 1} = idx - 1;
        day_length = unique(esc.hrs_map{idx,'hrs_per_day'});
        assert(size(day_length,1) == 1, 'Non-unique hours_per_day in day %s', num2str(day))
        mpc.days{day,2} = day_length;
    end
else
    mpc.days = [];
end

%%

