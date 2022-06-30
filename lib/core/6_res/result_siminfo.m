function sim_info = result_siminfo(esc, setup_opt)
% result_siminfo gathers information about simulation

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

% Sim Info
for field = fieldnames(esc)'
    f = char(field);
    c = class(esc.(f));
    if any(strcmp(c, {'double', 'char'}))
        if size(esc.(f), 1) == 1
            tmp.(f) = esc.(f);
        elseif size(esc.(f), 1) > 1 && size(esc.(f), 2) == 1
            tmp.(f) = esc.(f)';
        end
    end
end
sim_info.esc = struct2table(tmp);
clear tmp

for field = fieldnames(setup_opt)'
    f = char(field);
    c = class(setup_opt.(f));
    if any(strcmp(c, {'double', 'char'}))
        if size(setup_opt.(f), 1) == 1
            tmp.(f) = setup_opt.(f);
        elseif size(setup_opt.(f), 1) > 1 && size(setup_opt.(f), 2) == 1
            tmp.(f) = setup_opt.(f)';
        end
    end
end
sim_info.setup = struct2table(tmp);
clear tmp
