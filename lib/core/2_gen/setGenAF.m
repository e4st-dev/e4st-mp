function mpc = setGenAF(mpc, esc, eopt)
%setGenAF: Update generator availability factors (AF)

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Prepare Gen AF Data
info = esc.gen_af;
info = filterInfo(info, esc);

nh = size(esc.hrs_map, 1);

if isfield(info, 'joint')
    info_list = info.joint{:, :};
else
    info_list = (1:height(info.value))';
end

%% Initialize/Reset
% This section of the code creates the field mpc.availability_factor should
% it not already exist. Sets default availability factor to one. 
% This lists the availability factor for every generator in each hour
if isfield(mpc, 'availability_factor')
    new_gens = size(mpc.gen, 1) - size(mpc.availability_factor, 1);% a value
    mpc.availability_factor = [mpc.availability_factor; ones(new_gens, nh)];
    idx_update_af = mpc.newgen == max(mpc.newgen); % only update once
else
    mpc.availability_factor = ones(size(mpc.gen, 1), nh);
    idx_update_af = ones(size(mpc.gen, 1), 1) == 1;
end
mpc.idx_update_af = idx_update_af; % Save for tgtGenAF; %list of which gens to update


idx_updated = zeros(size(mpc.gen,1), 1); %create counter list to mark which generators were updated
%% Set Gen AF
for i = unique(info_list)'
    idx_info = i == info_list;
    cur_info = filterStruct(info, idx_info);
    info_tbl = info2table(cur_info);
    
    map_area = unique(cur_info.area{:, 'area'});
    assert(numel(map_area) == 1, 'Joint not properly upated in esc.gen_af. More than one area detected')
    gen_tbl = getMap(mpc, esc, 'gen', []);%list of all generators in gen and their areas
    gen_tbl = gen_tbl(:, {'bus', char(map_area), 'match_af'});
    gen_tbl{:, 'order'} = (1:height(gen_tbl))';

    %make class of generator map and availabilty factors areas the same
    if ~strcmp(class(info_tbl{:, 'subarea'}), class(gen_tbl{:, map_area}))
        if isa(gen_tbl{:, map_area}, 'double') && ~isa(info_tbl{:, 'subarea'}, 'double')
            gen_tbl = [gen_tbl, array2table(sprintfc('%d', gen_tbl{:, map_area}), 'VariableNames', {'subarea'})];
        end
    else
        gen_tbl{:, 'subarea'} = gen_tbl{:, map_area};
    end
    
    %AF map for every bus. 
    %Buses not in relevant area have NaN AF
    tmp_afs = outerjoin(gen_tbl, info_tbl, 'Type', 'Left', 'Keys', {'subarea', 'match_af'}, 'MergeKeys', true);
    %tmp_afs = tmp_afs(~isnan(tmp_afs{:,'bus'}),:);
    tmp_afs = sortrows(tmp_afs, 'order'); %availability factors

    if height(tmp_afs) ~= height(gen_tbl)
        error('Error in applying bus availability factors\n. More than one AF was assigned to a generator. ')
    end

    idx_gen = ~isnan(tmp_afs{:, end}); %keeps only generators that are at a substation covered by given joint in gen_af are at a substationo

    fuel = unique(info_tbl{:, 'genfuel'});
    type = unique(info_tbl{:, 'gentype'});
    if ~strcmp(fuel, '')
        idx_gen = idx_gen & strcmp(mpc.genfuel, fuel);
    end
    if ~strcmp(type, '')
        idx_gen = idx_gen & strcmp(mpc.gentype, type);
    end

    idx_gen = idx_gen & idx_update_af;

    if ~any(idx_gen)
        continue
    end
    
    %set availability factors
    s_idx = find(strcmp(tmp_afs.Properties.VariableNames, 'C1'));
    mpc.availability_factor(idx_gen, :) = tmp_afs{idx_gen, s_idx:(s_idx + nh - 1)};

    %mark as updated
    idx_updated = idx_updated | idx_gen;
    
    vfprintf(eopt.verbose, 'Availability factors set for %d generators by %s area and genfuel %s\n', sum(idx_gen), char(map_area), char(fuel));
end

%% check that all availability factors are plausible
assert(all(mpc.availability_factor >= 0, 'all'), "Generator availability factor less than zero")
assert(all(mpc.availability_factor <= 1, 'all'), "Generator availability factor greater than one")

%% check for any missing availability factors
%by default, availability factors for new units are set to one.
%This code thus looks for any generator where all availability factors in
%every hour are one.
% tol = 10^-5; 
% idx = abs(mean(mpc.availability_factor, 2)- 1) < tol;

%find indices of generators which were not updated
idx = ~idx_updated & mpc.idx_update_af;
idx = idx & ~strcmp(mpc.gentype, 'dl') & ~strcmp(mpc.gentype, 'test'); %make exception for test generators and loads.
idx = idx & mpc.gen(:, 8) == 1; %include only online generators

if sum(idx) > 0
    gentype = mpc.gentype{find(idx,1)};
    error(['A ' gentype ' unit has no assigned availability factor.'])
end    
    

% if strcmp(eopt.storage, 'T')
%     idx_dl = strcmp(mpc.genfuel, 'dl');
%     mpc.availability_factor(idx_dl, :) = mpc.availability_factor(idx_dl, :) .* esc.hrs_map{:, 'probability'}';
% end