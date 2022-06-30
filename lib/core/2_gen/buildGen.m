function [mpc, offer] = buildGen(mpc, esc, offer, eopt)
%buildGen: Create buildable generators for buildable generator type

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
% and Christoph Funke (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Pre-Process Data
info = esc.gen_build;
info = filterInfo(info, esc);

if isfield(info, 'joint')
    info_list = info.joint{:, :};
else
    info_list = (1:height(info.value))';
end

%% Process Data
define_constants;
for i = unique(info_list)'
    idx_info = i == info_list;
    cur_info = filterStruct(info, idx_info);
    info_tbl = info2table(cur_info);
    
    type = unique(info_tbl{:, 'gentype'});
    assert(length(type) == 1, 'More than one gentype contained in Joint number %u in gen build', i);
    % if the buildables are DAC, get treated differently below.
    is_dac = strcmp(unique(info_tbl{:, 'genfuel'}), 'dac');

    map_area = unique(cur_info.area{:, 'area'});
    %gen_tbl= getMap(mpc, esc, 'bus', map_area);
    gen_tbl = getMap(mpc, esc, 'bus', []);
    gen_tbl{:, 'order'} = (1:height(gen_tbl))';

    if ~strcmp(class(info_tbl{:, 'subarea'}), class(gen_tbl{:, map_area}))
        if isa(gen_tbl{:, map_area}, 'double') && ~isa(info_tbl{:, 'subarea'}, 'double')
            gen_tbl = [gen_tbl, array2table(sprintfc('%d', gen_tbl{:, map_area}), 'VariableNames', {'subarea'})];
        end
    else
        gen_tbl{:, 'subarea'} = gen_tbl{:, map_area};
    end

    %newgenstbl = outerjoin(gen_tbl, info_tbl, 'Type', 'Left', 'Keys', {'subarea'}, 'MergeKeys', true);
    newgenstbl = outerjoin(gen_tbl, info_tbl, 'Type', 'Full', 'Keys', {'subarea'}, 'MergeKeys', true);
    newgenstbl = newgenstbl(~isnan(newgenstbl{:, 'bus'}), :);
    newgenstbl = sortrows(newgenstbl, 'order');
    

    %if height(newgenstbl) ~= height(gen_tbl)
    %    error('Error in building new generators \n')
    %send

    idx_gen = ~isnan(newgenstbl{:, 'status'});
    

    if ~any(idx_gen)
        continue
    end

    newgenstbl = newgenstbl(idx_gen, :);

    numBus = sum(idx_gen);

    %% Gen
    newgens.gen = newgenstbl(:, 'bus');
    if is_dac
        newgens.gen{:, 'PG'} = -newgenstbl{:, 'CAPACITY'}/(1- esc.dist_loss); % dac implied distribution loss.
    else
        newgens.gen{:, 'PG'} = newgenstbl{:, 'CAPACITY'};
    end
    newgens.gen{:, 'QG'} = 0;
    newgens.gen{:, 'QMAX'} = newgens.gen{:, 'PG'} / 2;
    newgens.gen{:, 'QMIN'} = -newgens.gen{:, 'PG'} / 2;
    newgens.gen{:, 'VG'} = 1;
    newgens.gen{:, 'MBASE'} = 100;
    newgens.gen{:, 'GEN_STATUS'} = 1;
    if is_dac
        newgens.gen{:, 'QMAX'} = 0;
        newgens.gen{:, 'QMIN'} = 0;
        newgens.gen{:, 'PMAX'} = min(10, abs(newgens.gen{:, 'PG'}));
        newgens.gen{:, 'PMIN'} = newgens.gen{:, 'PG'};
    else
        newgens.gen{:, 'PMAX'} = newgens.gen{:, 'PG'};
        newgens.gen{:, 'PMIN'} = 0;
    end
    newgens.gen{:, 'PC1'} = 0;
    newgens.gen{:, 'PC2'} = 0;
    newgens.gen{:, 'QC1MIN'} = 0;
    newgens.gen{:, 'QC1MAX'} = 0;
    newgens.gen{:, 'QC2MIN'} = 0;
    newgens.gen{:, 'QC2MAX'} = 0;
    newgens.gen{:, 'RAMP_AGC'} = 0;
    newgens.gen{:, 'RAMP_10'} = Inf;
    newgens.gen{:, 'RAMP_30'} = 0;
    newgens.gen{:, 'RAMP_Q'} = 0;
    newgens.gen{:, 'APF'} = 0;
    newgens.gen = newgens.gen{:, :};

    %% Gencost
    newgens.gencost = array2table(repmat(2, numBus, 1), 'VariableNames', {'MODEL'});
    newgens.gencost{:, 'STARTUP'} = 0;
    newgens.gencost{:, 'SHUTDOWN'} = 0;
    newgens.gencost{:, 'NCOST'} = 2;
    newgens.gencost{:, 'GENCOST'} = newgenstbl{:, 'VOM'}; % Gets updated properly later in updateGenCost with fuel cost (and DAC updated with proper VOM)
    newgens.gencost = newgens.gencost{:, :};
    newgens.gencost = [newgens.gencost, zeros(numBus, size(mpc.gencost, 2) - 5)];

    %% Genfuel/Gentype/Newgen
    newgens.genfuel = newgenstbl{:, 'genfuel'};
    newgens.gentype = newgenstbl{:, 'gentype'};
    newgens.newgen = ones(numBus, 1);

    %% Get map
    %loc_map = getMap(mpc, esc, 'bus', []); loc_map = loc_map(idx_bus, :);

    %% Gen_aux
    newgens.gen_aux = newgenstbl(:, 'HR');
    if is_dac
        % HR controls how much NG bought per MWh dac usage
        newgens.gen_aux{:, 'HR'} = - newgens.gen_aux{:, 'HR'}*(1-esc.dist_loss);
        newgens.gen_aux{:, 'EMIS_CO2'} = newgenstbl{:, 'EMIS_CO2'}*(1-esc.dist_loss);
        newgens.gen_aux{:, 'EMIS_NOX'} = newgenstbl{:, 'EMIS_NOX'}*(1-esc.dist_loss);
        newgens.gen_aux{:, 'EMIS_SO2'} = newgenstbl{:, 'EMIS_SO2'}*(1-esc.dist_loss);
        newgens.gen_aux{:, 'EMIS_PM25'} = newgenstbl{:, 'EMIS_PM25'}*(1-esc.dist_loss);
    else
        newgens.gen_aux{:, 'EMIS_CO2'} = newgenstbl{:, 'EMIS_CO2'};
        newgens.gen_aux{:, 'EMIS_NOX'} = newgenstbl{:, 'EMIS_NOX'};
        newgens.gen_aux{:, 'EMIS_SO2'} = newgenstbl{:, 'EMIS_SO2'};
        newgens.gen_aux{:, 'EMIS_PM25'} = newgenstbl{:, 'EMIS_PM25'};
    end
    newgens.gen_aux{:, 'EMIS_CO2e'} = newgenstbl{:, 'EMIS_CO2e'}; % gets written over with calculation of CO2e in updateCO2e
    newgens.gen_aux{:, 'DEATHS_INF_NOX'} = 0; %newgenstbl{:, 'DEATHS_INF_NOX'};
    newgens.gen_aux{:, 'DEATHS_ADLT_NOX'} = 0; %newgenstbl{:, 'DEATHS_ADLT_NOX'};
    newgens.gen_aux{:, 'DEATHS_INF_SO2'} = 0; %newgenstbl{:, 'DEATHS_INF_SO2'};
    newgens.gen_aux{:, 'DEATHS_ADLT_SO2'} = 0; %newgenstbl{:, 'DEATHS_ADLT_SO2'};
    newgens.gen_aux{:, 'VSL_INF'} = 0;
    newgens.gen_aux{:, 'VSL_ADLT'} = 0;
    newgens.gen_aux{:, 'DEATHS_NOX_US_LO'} = 0;
    newgens.gen_aux{:, 'DEATHS_SO2_US_LO'} = 0;
    newgens.gen_aux{:, 'DEATHS_NOX_US_HI'} = 0;
    newgens.gen_aux{:, 'DEATHS_SO2_US_HI'} = 0;
    newgens.gen_aux{:, 'DAM_CO2'} = 0;
    newgens.gen_aux{:, 'DAM_CH4'} = 0;
    newgens.gen_aux{:, 'ESH_TYPE'} = newgenstbl{:, 'ESH_TYPE'};
    newgens.gen_aux{:, 'CH4_FUEL_CONTENT'} = 0;
    newgens.gen_aux{:, 'FUEL_COST'} = 0;
    if is_dac %variable costs for DAC have to be negative, since power is negative
        newgens.gen_aux{:, 'VOM'} = -newgenstbl{:, 'VOM'}*(1-esc.dist_loss);
        newgens.gen_aux{:, 'FOM'} = newgenstbl{:, 'FOM'}*(1-esc.dist_loss);
        newgens.gen_aux{:, 'CAP_COST'} = newgenstbl{:, 'CAP_COST'}*(1-esc.dist_loss);
    else
        newgens.gen_aux{:, 'VOM'} = newgenstbl{:, 'VOM'};
        newgens.gen_aux{:, 'FOM'} = newgenstbl{:, 'FOM'};
        newgens.gen_aux{:, 'CAP_COST'} = newgenstbl{:, 'CAP_COST'};
    end
    newgens.gen_aux{:, 'TRANS_COST'} = newgenstbl{:, 'TRANS_COST'};
    newgens.gen_aux{:, 'ROUTINE_CAPEX'} = 0;
    newgens.gen_aux{:, 'PTC'} = 0;
    newgens.gen_aux{:, 'ITC'} = 0;
    newgens.gen_aux{:, 'PAST_CAPEX'} = newgenstbl{:, 'CAP_COST'} + newgenstbl{:, 'TRANS_COST'}; 
    newgens.gen_aux{:, 'CHP_CO2_MULTI'} = 1;
    newgens.gen_aux{:, 'CAP_CREDIT'} = 0;
    newgens.gen_aux{:, 'HIST_CF'} = 0;
    newgens.gen_aux{:, 'TGT_AF'} = 0;
    newgens.gen_aux{:, 'TGT_CF'} = 0;
    newgens.gen_aux{:, 'CO2_CAPTURE_PC'} = newgenstbl{:, 'CO2_CAPTURE_PC'};
    % Use capture percentage and emission rate (after capture) to calculate how much is stored
    % This is not built for capture Pctage of 1 or negative emissions.
    newgens.gen_aux{:, 'STORAGE_CO2'} = newgens.gen_aux{:, 'EMIS_CO2'} ./ (1 - newgens.gen_aux{:, 'CO2_CAPTURE_PC'}) .* newgens.gen_aux{:, 'CO2_CAPTURE_PC'};
    if is_dac
        % DAC "generation" is negative, so multiplying negative power by
        % negative storage will yield positive storage.
        % Capture_pc has a special meaning for dac, it is additional CO2
        % produced and immediately captured in process.
        newgens.gen_aux{:, 'CO2_CAPTURE_PC'} = newgenstbl{:, 'CO2_CAPTURE_PC'}*(1-esc.dist_loss);
        newgens.gen_aux{:, 'STORAGE_CO2'} = -(newgens.gen_aux{:, 'EMIS_CO2'} + newgens.gen_aux{:, 'CO2_CAPTURE_PC'});
    end
    
    newgens.gen_aux{:, 'STORAGE_EFFICIENCY'} = newgenstbl{:, 'STORAGE_EFFICIENCY'};
    newgens.gen_aux{:, 'STORAGE_E_CAPACITY'} = newgenstbl{:, 'STORAGE_E_CAPACITY'};
    
    %% Storage
    nge = size(mpc.gen,1); %number of existing generators
    ngn = size(newgenstbl,1);%number of new generators
    newgens.short_term_storage = [(nge+1):(nge+ngn)]'; %new generator indices
    newgens.short_term_storage = ...
        [newgens.short_term_storage, newgenstbl{:, 'STORAGE_EFFICIENCY'}];
    newgens.short_term_storage = ...
        [newgens.short_term_storage, newgenstbl{:, 'STORAGE_E_CAPACITY'}];
    

    %% Gen_desc
    if isfield(cur_info, 'type')
        newgens.gen_desc = array2table(repmat({strjoin([cur_info.gentype{:, :} cur_info.type{:, :}], {'_'})}, numBus, 1), 'VariableNames', {'NAME'});
    else
        newgens.gen_desc = array2table(repmat(unique(cur_info.gentype{:, :}), numBus, 1), 'VariableNames', {'NAME'});
    end
    newgens.gen_desc{:, 'PLANTID'} = {'na'};
    newgens.gen_desc{:, 'GENID'} = {'na'};
    newgens.gen_desc{:, 'ON_YR'} = esc.year;
    newgens.gen_desc{:, 'ON_SIMYR'} = {eopt.sim_yr};
    newgens.gen_desc{:, 'RET_YR'} = 9999;
    newgens.gen_desc{:, 'RET_SIMYR'} = {'na'};
    newgens.gen_desc{:, 'P_RET_YR'} = 9999;
    newgens.gen_desc{:, 'RET_REASON'} = {'na'};
    newgens.gen_desc{:, 'INVEST_TYPE'} = newgenstbl{:, 'INVEST_TYPE'};
    newgens.gen_desc{:, 'AVG_CAP'} = newgenstbl{:, 'CAPACITY'};
    newgens.gen_desc{:, 'HIST_CAP'} = 0;
    newgens.gen_desc{:, 'HIST_GEN'} = 0;
    newgens.gen_desc{:, 'HIST_CF'} = 0;
    newgens.gen_desc{:, 'AGE'} = 0;
    newgens.gen_desc{:, 'CHP'} = {'No'};

    %% Gen_map
    count = 0;
    for area = mpc.gen_map.Properties.VariableNames
        if ismember(area, {'gen_id', 'reg_factor_prelim', 'reg_factor'})
            continue
        elseif ~ismember(area, newgenstbl.Properties.VariableNames)
            newgenstbl{:, area} = {'err_tmp'}; %set default if not found
        end
        count = count + 1;
        if count == 1
            newgens.gen_map = newgenstbl(:, area);
        else
            newgens.gen_map{:, area} = newgenstbl{:, area};
        end
    end
    newgens.gen_map{:, 'reg_factor'} = newgens.gen_map{:, 'state_reg_multi'};
    newgens.gen_map{:, 'reg_factor_prelim'} = 0; %prelim is just for pre-existing units
    newgens.gen_map{:, 'gen_id'} = ((max(mpc.gen_map{:, 'gen_id'}) + 1):(max(mpc.gen_map{:, 'gen_id'}) + numBus))';

    %% Offer
    newgens.offer = zeros(numBus, 13);
    if is_dac
        newgens.offer(:, 4) = abs(newgens.gen(:, PMIN)); % Installed Cap
    else
        newgens.offer(:, 2) = newgens.gen(:, PMAX); % Installed Cap
        newgens.offer(:, 4) = Inf;
    end


    %% Add to MPC and Offer
    [mpc, offer] = append_e4st_gen(mpc, offer, newgens);

    vfprintf(eopt.verbose, '%d buildable %s generators totalling %.1f capacity are added by %s area\n', numBus, char(unique(cur_info.gentype{:, :})), sum(newgenstbl{:, 'CAPACITY'}), char(map_area));
end

mpc.ng = size(mpc.gen, 1);
mpc.nb = size(mpc.bus, 1);

idx_gen = mpc.newgen == 1;
offer = updateOfferPrc(mpc, offer, idx_gen, eopt);
mpc = updateGenCost(mpc, idx_gen, eopt);
