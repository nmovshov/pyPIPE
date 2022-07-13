function options = pipeset(varargin)
%PIPESET Create options structure used by PIPE classes.
%   OPTIONS = PIPESET('NAME1',VALUE1,'NAME2',VALUE2,...) creates an options
%   structure OPTIONS in which the named properties have the specified values. Any
%   unspecified properties have default values. Case is ignored for property
%   names.
%
%   PIPESET with no input arguments displays all property names and their
%   possible values.
%
%KNOWN PROPERTIES
%
%reducelinks - Convert from *Planet to struct on the fly to save memory [ {0} 1 2 3 ]
%backupfreq - Backup frequency (save temp file every n links) [ nonnegative integer {0} ]
%verbosity - Level of runtime messages [0 {1} 2 3]
%linkJ0 - Use Js of previous link to speed up convergence of new link [ true | {false} ]
%quickreject - Attempt to quick-reject model without converging [ {true} | false ]
%matchmass - Modify converged models' density to match observables.M [ true | {false} ]
%matchreq - Modify converged models' level radii to match observables.a0 [ true | {false} ]

% If no arguments print usage and return.
if (nargin == 0) && (nargout == 0)
    print_usage()
    return
end

% Define name-value pairs.
p = inputParser;
p.FunctionName = mfilename;

p.addParameter('reducelinks',0,@isnonnegintscalar)
p.addParameter('linkJ0',false,@islogicalscalar)
p.addParameter('quickreject',true,@islogicalscalar)
p.addParameter('backupfreq',0,@isnonnegintscalar)
p.addParameter('verbosity',1,@isnonnegintscalar)
p.addParameter('matchmass',false,@islogicalscalar)
p.addParameter('matchreq',false,@islogicalscalar)

% Parse name-value pairs and return.
p.parse(varargin{:})
options = p.Results;

end

function islogicalscalar(x)
validateattributes(x,{'logical'},{'scalar'})
end

function isfunctionhandle(x)
validateattributes(x,{'function_handle'},{})
end

function isnonnegintscalar(x)
validateattributes(x,{'numeric'},{'nonnegative','integer','scalar'})
end

function print_usage()
help(mfilename)
end
