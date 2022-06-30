function info_tbl = info2table(info)
% info2table Convert esc info struct into a table

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Convert esc info struct into a table
info_tbl = [];
for field = fieldnames(info)'
	info_tbl = [info_tbl info.(char(field))];
end
