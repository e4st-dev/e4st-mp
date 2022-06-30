function rv = e4st_ver(varargin)
%E4ST_VER  Prints or returns E4ST version info for current installation.
%   V = E4ST_VER returns the current E4ST version number.
%   V = E4ST_VER('all') returns a struct with the fields Name, Version,
%   Release and Date (all strings). Calling E4ST_VER without assigning the
%   return value prints the version and release date of the current
%   installation of E4ST.
%
%   See also MPVER.

%   E4ST
%   Copyright (c) 2010-2022 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).

v = struct( 'Name',     'E4ST', ...
            'Version',  '2.1b1', ...
            'Release',  '', ...
            'Date',     '23-Jul-2020');
if nargout > 0
    if nargin > 0
        rv = v;
    else
        rv = v.Version;
    end
else
    fprintf('%-22s Version %-9s  %11s\n', v.Name, v.Version, v.Date);
end
