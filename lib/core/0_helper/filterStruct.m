function struct = filterStruct(struct, idx_keep, fields)
% filterStruct Filter rows of fields in a struct
% - idx_keep: specifies index of rows to keep
% - fields: specifies fields of struct to filter. if not provided, all
%           fields are filtered

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% All fields if not specified
if nargin < 3
    fields = fieldnames(struct)';
end

%% Filter struct
for field = fields
    try
    %if size(struct.(char(field)), 1) == size(idx_keep, 1)
        struct.(char(field)) = struct.(char(field))(idx_keep, :);
    %else
    catch
        disp('Warning: idx_keep not same size as struct');
    end
end
