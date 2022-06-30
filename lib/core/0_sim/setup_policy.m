function [mpc, offer] = setup_policy(mpc, offer, esc, eopt)
% SETUP_POLICY applies power sector policies to a simulation

% E4ST 
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)

% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

vfprintf(eopt.verbose, '\n\n-----Running E4ST Set-Up Module: Policies-----\n\n')

%% Initalize
% set policy spreadsheet to info
info = esc.policy;

% the filterInfo function filters to  policy file to only include policies where
% status = 1, the start_year is before simulation year, end year is after
% simulation, and filters relevant year
info = filterInfo(info, esc);

%filter value to only include year in question
info.value = info.value(:, ['Y', num2str(esc.year)]);

%info.joint assigns a unique number to each policy
if isfield(info, 'joint')
    info_list = info.joint{:, :};
else
    info_list = (1:height(info.value))';
end

%% Load (For RPS, RPS_QTY, IRM)
%if any(strcmp(info.pol{:, 'type'}, 'rps') | strcmp(info.pol{:, 'type'}, 'rps_qty') | strcmp(info.pol{:, 'type'}, 'irm'))
% scales load to match growth rates and such. load info needed for some policies
[tmp_mpc, ~] = setup_load(mpc, esc, eopt, 0);

% Things which have to pay for credits for elec they consume (demand side)
idx_ps_dl = ismember(info.genfuel{:, :}, {'dl', 'dac', 'storage'}) & ismember(info.pol{:, 'type'}, ...
    {'rps', 'ces', 'ces_coal', 'ces_ngcc'});

% For elec consumers, this section brings the pctage they have to pay into gen_wgt
info.gen_wgt{idx_ps_dl, :} = info.value{idx_ps_dl, :};

%% Apply Policies

%for all individual policies ....
for i = unique(info_list)'
    %identify rows corresponding to that policy :)
    idx_info = i == info_list;

    %filter results to rows relevant for the policy
    cur_info = filterStruct(info, idx_info);

    %specify which generators are affected
    %creates index of generators in mpc.gen that satisifes all criteria of
    %the policy, specifically: area,
    %note idx_gen and gen_wgt are the same, unless info.gen_wgt is not
    %equal to one for some generator type
    [idx_gen, gen_wgt] = getInfoIdx(mpc, esc, cur_info, 'gen');

    if ~any(idx_gen)
        continue
    end

    pol_type = char(unique(cur_info.pol{:, 'type'}));
    pol_name = char(unique(cur_info.pol{:, 'name'}));
    pol_val = unique(cur_info.value{:, :}); %policy value in given year

    %go through the different policy types and apply the relevant changes
    switch pol_type
        case 'calib_qty'
            mpc = applyCalib(mpc, esc, pol_val, 'qty', idx_gen, pol_name, eopt);
        case 'calib_prc'
            mpc = applyCalib(mpc, esc, pol_val, 'prc', idx_gen, pol_name, eopt);
        case 'ptc'
            mpc = applyPTC(mpc, pol_val, idx_gen, eopt);
        case 'itc'
            [mpc, offer] = applyITC(mpc, offer, pol_val, idx_gen, eopt);
        case 'irm'
            assert(~strcmp(eopt.dac, 'T'), "IRM policies not set up for use with DAC yet");
            if pol_val <= 1 % percentage, convert to level
                idx_bus = getInfoIdx(mpc, esc, cur_info, 'bus');
                pol_val = max(sum(tmp_mpc.hourly_load(idx_bus, :)), [], 2) * (1 + pol_val);
            end
            idx_gen = ones(length(idx_gen), 1) .* idx_gen .* mpc.gen_aux{:, 'CAP_CREDIT'};
            mpc = applyIRMReq(mpc, pol_val, idx_gen, pol_name, eopt);
        case 'uplift'
            assert(~strcmp(eopt.dac, 'T'), "Uplift policies not set up for use with DAC yet");
            mpc = applyUplift(mpc, pol_val, idx_gen);
        case 'emis_prc_co2'
            mpc = applyEmisPol(mpc, pol_val, 'CO2', 'prc', idx_gen, pol_name, esc, eopt);
        case 'emis_prc_co2e'
            mpc = applyEmisPol(mpc, pol_val, 'CO2e', 'prc', idx_gen, pol_name, esc, eopt);
        case 'emis_prc_co2e_20y'
            mpc = applyEmisPol(mpc, pol_val, 'CO2e_20y', 'prc', idx_gen, pol_name, esc, eopt);
        case 'emis_prc_nox'
            mpc = applyEmisPol(mpc, pol_val, 'NOX', 'prc', idx_gen, pol_name, esc, eopt);
        case 'emis_prc_so2'
            mpc = applyEmisPol(mpc, pol_val, 'SO2', 'prc', idx_gen, pol_name, esc, eopt);
        case 'emis_prc_pm'
            mpc = applyEmisPol(mpc, pol_val, 'PM25', 'prc', idx_gen, pol_name, esc, eopt);
        case 'emis_cap_co2'
            mpc = applyEmisPol(mpc, pol_val, 'CO2', 'cap', idx_gen, pol_name, esc, eopt);
        case 'emis_cap_co2e'
            mpc = applyEmisPol(mpc, pol_val, 'CO2e', 'cap', idx_gen, pol_name, esc, eopt);
        case 'emis_cap_nox'
            mpc = applyEmisPol(mpc, pol_val, 'NOX', 'cap', idx_gen, pol_name, esc, eopt);
        case 'emis_cap_so2'
            mpc = applyEmisPol(mpc, pol_val, 'SO2', 'cap', idx_gen, pol_name, esc, eopt);
        case 'co2_storage_cap'
            mpc = applyEmisPol(mpc, pol_val, 'CO2_STOR', 'cap', idx_gen, pol_name, esc, eopt);
        case 'co2_storage_prc'
            mpc = applyEmisPol(mpc, pol_val, 'CO2_STOR', 'prc', idx_gen, pol_name, esc, eopt);
        case 'rps'
            %returns index of relevant generators (idx_gen)
            %returns list with policy_value at generators of type 'dl'(idx_dl)
            %changes policy_value to 1

            [idx_gen, idx_dl, pol_val] = setupPS(tmp_mpc, esc, idx_gen, gen_wgt, cur_info, pol_val);
            %cut off policy if not properly specified
            if sum(idx_dl) == 0 || isnan(sum(pol_val))
                continue
            end
            %applies policy to mpc file
            mpc = applyPS(mpc, esc, pol_val, 'perc', 'rps', idx_gen, idx_dl, pol_name, eopt);
        case 'rps_qty'
            %idx_area_gen = getInfoIdx(mpc, esc, rmfield(cur_info, {'genfuel', 'gentype'}), 'gen');
            if pol_val >= 0 && pol_val <= 1 % target value, convert to total output quantity
                idx_area_bus = getInfoIdx(mpc, esc, rmfield(cur_info, {'genfuel', 'gentype'}), 'bus');
                pol_val = sum(tmp_mpc.ann_load(idx_area_bus)) * pol_val;
            end
            mpc = applyPS(mpc, esc, pol_val, 'qty', 'rps', idx_gen, [], pol_name, eopt);
        case 'rps_qty_max'
            %idx_area_gen = getInfoIdx(mpc, esc, rmfield(cur_info, {'genfuel', 'gentype'}), 'gen');
            if pol_val >= 0 && pol_val <= 1 % target value, convert to total output quantity
                idx_area_bus = getInfoIdx(mpc, esc, rmfield(cur_info, {'genfuel', 'gentype'}), 'bus');
                pol_val = sum(tmp_mpc.ann_load(idx_area_bus)) * pol_val;
            end
            mpc = applyPS(mpc, esc, pol_val, 'qty_max', 'rps', idx_gen, [], pol_name, eopt);
        case 'rps_prc'
            mpc = applyPS(mpc, esc, pol_val, 'prc', 'rps', idx_gen, [], pol_name, eopt);
        case 'ces'
            [idx_gen, idx_dl, pol_val] = setupPS(tmp_mpc, esc, idx_gen, gen_wgt, cur_info, pol_val);
            if sum(idx_dl) == 0 || isnan(sum(pol_val))
                continue
            end
            mpc = applyPS(mpc, esc, pol_val, 'perc', pol_name, idx_gen, idx_dl, pol_name, eopt);
        case 'ces_coal'
            [idx_gen, idx_dl, pol_val] = setupPS(tmp_mpc, esc, idx_gen, gen_wgt, cur_info, pol_val);
            if sum(idx_dl) == 0 || isnan(sum(pol_val))
                continue
            end
            mpc = applyPS(mpc, esc, pol_val, 'perc', 'ces_coal', idx_gen, idx_dl, pol_name, eopt);
        case 'ces_ngcc'
            [idx_gen, idx_dl, pol_val] = setupPS(tmp_mpc, esc, idx_gen, gen_wgt, cur_info, pol_val);
            if sum(idx_dl) == 0 || isnan(sum(pol_val))
                continue
            end
            mpc = applyPS(mpc, esc, pol_val, 'perc', 'ces_ngcc', idx_gen, idx_dl, pol_name, eopt);
        case 'ces_qty'
            if pol_val <= 1 % target value, convert to total output quantity
                idx_area_bus = getInfoIdx(mpc, esc, rmfield(cur_info, {'genfuel', 'gentype'}), 'bus');
                pol_val = sum(tmp_mpc.ann_load(idx_area_bus)) * pol_val;
            end
            mpc = applyPS(mpc, esc, pol_val, 'qty', 'ces', idx_gen, [], pol_name, eopt);
        case 'ces_prc'
            mpc = applyPS(mpc, esc, pol_val, 'prc', 'ces', idx_gen, [], pol_name, eopt);
        case 'ces_prc_ngcc'
            mpc = applyPS(mpc, esc, pol_val, 'prc', 'ces_ngcc', idx_gen, [], pol_name, eopt);            
        otherwise
            vfprintf(eopt.verbose, 'Policy type %s not recognized\n', char(pol_type))
    end
end

%% Past Capex
% For buildable generators, PAST_CAPEX stores invest cost. Apply ITC here.
idx_new = mpc.newgen == 1;
mpc.gen_aux{idx_new, 'PAST_CAPEX'} = mpc.gen_aux{idx_new, 'CAP_COST'} + mpc.gen_aux{idx_new, 'TRANS_COST'} - mpc.gen_aux{idx_new, 'ITC'};


end
