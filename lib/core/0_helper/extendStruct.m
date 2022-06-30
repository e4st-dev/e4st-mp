function [struct] = extendStruct(struct, struct_add, fields)
% extendStruct Stacks two different data structures with the same field names.
% The entries in 'struct_add' are appended at the end of
% the corresponding fields in 'struct'
%INPUTS:
%  struct ... the data structure to be appended to
%  struct_add ... the data structure containing the values to be appended to the fields in struct
%  fields (optional) ... list of fields to be extended. If none is
%  provided, then all fields are used by default.

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

if nargin < 3
    fields = fieldnames(struct);
    idx = ismember(fieldnames(struct), fieldnames(struct_add));
    fields = fields(idx,:);
end

for i = 1:length(fields)
    field = fields{i};
    struct.(field) = [struct.(field); struct_add.(field)];
end
end

