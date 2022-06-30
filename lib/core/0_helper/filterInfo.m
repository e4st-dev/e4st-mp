function info = filterInfo(info, esc)
% filterInfo Selects active E4ST setup elements
% Select if...
% - STATUS: Equals 1
% - START_YEAR: Less than or equal to simulation year
% - END_YEAR: Greater than or equal to simulation year
% - YEAR: Equals simulation year

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Filter info
% If status equals 1
if isfield(info, 'status')
	idx_info = info.status{:, :} == 1;
	info = filterStruct(info, idx_info);
end

% If start_year less than or equal to simulation year
if isfield(info, 'start_year')
	idx_info = esc.year >= info.start_year{:, :};
	info = filterStruct(info, idx_info);
end

% If end_year less than or equal to simulation year
if isfield(info, 'end_year')
	idx_info = esc.year <= info.end_year{:, :};
	info = filterStruct(info, idx_info);
end

% If year equals simulation year
if isfield(info, 'year')
	idx_info = esc.year == info.year{:, :};
	info = filterStruct(info, idx_info);
end
