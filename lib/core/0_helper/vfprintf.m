function vfprintf(verbose, varargin)
% vfprintf Display progress information if verbose = 1

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Display progress information if verbose = 1
if isstruct(verbose)
    %save logs to file
    file_id = fopen(verbose.filename, 'a');
    fprintf(file_id, varargin{:});
    fprintf(file_id, '\n');
    fclose(file_id);
    
    %print to command window
    if verbose.verbose
        fprintf(varargin{:})     
    end    
else
    if verbose
        fprintf(varargin{:})
    end    
end