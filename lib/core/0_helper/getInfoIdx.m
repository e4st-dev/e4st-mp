function [idx_all, wgt_all] = getInfoIdx(mpc, esc, info, type, postRet)
% getInfoIdx Return specified index for buses, generators, or branches
% - idx_all: index of 0 or 1
% - wgt_all: weight between 0 and 1
%
% For type = "gen", returns the indices of the generators that are included
% in the policy. Also uses the gen_wgt attribute of the policy to set
% weights for the generators in wgt_all.
% Optional postRet = 1 will include gens that were endogenously retired in
% the current year. Meant to be used in postprocessing.

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Argument Parsing
if nargin < 5
    postRet = 0;
end

%% Initialize index
define_constants;
if strcmp(type, 'gen')
    idx_all = zeros(size(mpc.gen, 1), 1);
elseif strcmp(type, 'bus')
    idx_all = zeros(size(mpc.bus, 1), 1);
elseif strcmp(type, 'branch')
    idx_all = zeros(size(mpc.branch, 1), 1);
end
wgt_all = ones(size(idx_all, 1), 1);

%% Get index
for i = 1:height(info.value)

    idx = ones(size(idx_all, 1), 1);
    wgt = ones(size(idx_all, 1), 1);

    for opt = fieldnames(info)'
        if ~strcmp(type, 'gen') && ~strcmp(opt, 'area')
            continue
        end

        switch char(opt)
            case 'area'
                tmp_area = info.area{i, 'area'};
                if ~strcmp(tmp_area, '')
                    if strcmp(type, 'branch')
                        mpc.bus = mpc.branch;
                        loc_map = getMap(mpc, esc, 'bus', tmp_area);
                    else
                        loc_map = getMap(mpc, esc, type, tmp_area);
                    end

                    tmp_map = loc_map{:, tmp_area};
                    tmp_loc = info.area{i, 'subarea'};
                    if isa(tmp_map, 'double') && isa(tmp_loc, 'double')
                        idx = idx & tmp_map == tmp_loc;
                    elseif isa(tmp_map, 'double')
                        %tmp_map = num2str(tmp_map);
                        tmp_map = sprintfc('%d',tmp_map);
                        idx = idx & strcmp(tmp_map, tmp_loc);
                    elseif isa(tmp_loc, 'double')
                        tmp_loc = num2str(tmp_loc);
                        idx = idx & strcmp(tmp_map, tmp_loc);
                    else
                        idx = idx & strcmp(tmp_map, tmp_loc);
                    end
                end
            case 'genfuel'
                if ~strcmp(info.genfuel{i, :}, '') && strcmp(type, 'gen')
                    idx = idx & strcmp(mpc.genfuel, info.genfuel{i, :});
                end
            case 'gentype'
                if ~strcmp(info.gentype{i, :}, '') && strcmp(type, 'gen')
                    idx = idx & strcmp(mpc.gentype, info.gentype{i, :});
                end
            case 'newgen'
                if ~isnan(info.newgen{i, :}) && strcmp(type, 'gen')
                    idx = idx & mpc.newgen == 1;
                end
            case 'min_age'
                if ~isnan(info.min_age{i, :}) && strcmp(type, 'gen')
                    idx = idx & mpc.gen_desc{:, 'AGE'} >= info.min_age{i, :};
                end
            case 'max_age'
                if ~isnan(info.max_age{i, :}) && strcmp(type, 'gen')
                    idx = idx & mpc.gen_desc{:, 'AGE'} < (info.max_age{i, :} + esc.year_delta);
                end
            case 'min_cap'
                if ~isnan(info.min_cap{i, :}) && strcmp(type, 'gen')
                    idx = idx & mpc.gen_desc{:, 'AVG_CAP'} >= info.min_cap{i, :};
                end
            case 'max_cap'
                if ~isnan(info.max_cap{i, :}) && strcmp(type, 'gen')
                    idx = idx & mpc.gen_desc{:, 'AVG_CAP'} < info.max_cap{i, :};
                end
            case 'min_on_yr'
                if ~isnan(info.min_on_yr{i, :}) && strcmp(type, 'gen')
                    idx = idx & mpc.gen_desc{:, 'ON_YR'} >= info.min_on_yr{i, :};
                end
            case 'max_on_yr'
                if ~isnan(info.max_on_yr{i, :}) && strcmp(type, 'gen')
                    idx = idx & mpc.gen_desc{:, 'ON_YR'} < info.max_on_yr{i, :};
                end
        end
    end

    idx_all = idx_all | idx;
    
    if isfield(info, 'gen_wgt') && info.gen_wgt{i, :} ~= 1 && ~isnan(info.gen_wgt{i, :}) && strcmp(type, 'gen')
        wgt_all(idx) = wgt(idx) .* info.gen_wgt{i, :};
    end

end

% Only keep index for generators online
% If run after retireGen, want to include recently retired gens in results
if strcmp(type, 'gen')
    if postRet
        idx_ret = mpc.gen_desc{:, 'RET_YR'} == esc.year & strcmp(mpc.gen_desc{:,'RET_REASON'}, 'endogenous') ;
        idx_all = idx_all & (mpc.gen(:, GEN_STATUS) == 1 | idx_ret);
        wgt_all = wgt_all .* (idx_all == 1);        
    else
        idx_all = idx_all & mpc.gen(:, GEN_STATUS) == 1;
        wgt_all = wgt_all .* (idx_all == 1);
    end
end
