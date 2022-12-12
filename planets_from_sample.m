function planets = planets_from_sample(samplefile,observables,varargin)
%PLANETS_FROM_SAMPLE Convert a ppwd-sample to array of TOFPlanets
%   planets_from_sample(samplefile, observables, key1=value1, key2=value2,...)
%   reads sample-space values from samplefile and creates TOFPlanet objects
%   cooked to equilibrium. The array of TOFPlanets is saved to .mat file.
%
%KNOWN PARAMETER KEYS (case ignored and unique partial match allowed)
%
%nsamples - Number of planets to create, 0 for full sample [ positive integer {0} ]
%fixrot - If true use observables.P else use first element of sample-space vector [ {true} | false ]
%toforder - Theory of figures expansion order [ {4} 7 ]
%toflevels - Number of level surfaces used to discretize density [ positive integer {4096} ]
%xlevels - Skip-and-spline levels [ integer (-1 to disable) {256} ]
%savess - Retain full shape information (triples file size) [ true | {false} ]
%rdclvl - Reduction level [ 0=none, {1=to struct}, 2=to single, 3=to scalars ]
%savemat If true save to mat file [ true | {false} ]

% If no arguments print usage and return.
if (nargin == 0) && (nargout == 0)
    help(mfilename)
    return
end

% Inputs
narginchk(2,inf)
validateattributes(samplefile,{'char'},{})
validateattributes(observables,{'observables'},{})
args = PCL(varargin{:});

% Load sample-space array from file
sample = readmatrix(samplefile);

% Create a planet from each row
nsamp = args.nsamples;
if nsamp == 0, nsamp = size(sample,1); end
fprintf('Cooking planets...\n')
planets = [];
t0 = tic;
for k=1:nsamp
    fprintf('Cooking planet %d of %d...',k,nsamp)
    s = sample(k,:);
    p = cook_planet(s,observables,args);
    if args.rdclvl > 0
        p = p.to_struct(args.rdclvl,args.savess);
    end
    planets = [planets, p]; %#ok<AGROW> 
    fprintf('done.\n')
end
fprintf('Cooking planets...done. (%0.2g sec.)\n',toc(t0))

% Save cooked planets
if args.savemat
    outname = [samplefile(1:end-4), '_planets.mat'];
    save(outname,"planets")
    fprintf('mPickled %d planets in %s.\n',length(planets),outname)
end
end

%% The cooker
function tp = cook_planet(x,obs,args)
if args.fixrot
    Prot = obs.P;
else
    Prot = x(1)*obs.dP/2 + obs.P;
    x = x(2:end);
end
tp = ppwd.ppwd_planet(args.toflevels,x,obs);
tp.period = Prot;
tp.opts.toforder = args.toforder;
tp.opts.xlevels = args.xlevels;
tp.relax_to_rotation();
end

%% Helpers
function args = PCL(varargin)
p = inputParser;
p.FunctionName = mfilename;
p.addParameter('nsamples', 0, @isnonnegintscalar)
p.addParameter('fixrot',true, @islogicalscalar)
p.addParameter('toforder', 4, @(n)(n==4)||(n==7))
p.addParameter('toflevels',4096,@isposintscalar)
p.addParameter('xlevels',256,@isintscalar)
p.addParameter('savess',false,@islogicalscalar)
p.addParameter('rdclvl',1,@isnonnegintscalar)
p.addParameter('savemat',false,@islogicalscalar)

% Parse name-value pairs and return.
p.parse(varargin{:})
args = p.Results;
end

function isintscalar(x)
validateattributes(x,{'numeric'},{'integer','scalar'})
end

function isposintscalar(x)
validateattributes(x,{'numeric'},{'positive','integer','scalar'})
end

function isnonnegintscalar(x)
validateattributes(x,{'numeric'},{'nonnegative','integer','scalar'})
end

function islogicalscalar(x)
validateattributes(x,{'logical'},{'scalar'})
end
