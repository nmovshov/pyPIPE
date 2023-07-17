function savehc(fh, varargin)
%SAVEHC Save and export figure in triplocates.
%   SAVEHC(fh) saves the figure identified by handle fh to the file newfig.fig and
%   also exports it as newfig.png and newfig.pdf, all in the current directory.
%   The exports are done with export_fig, with -r600 and -nocrop.
%
%   SAVEHC(fh, 'name', value) recognizes optional parameters (partial match okay):
%       'name': string or char, base name of saved files (default: 'newfig')
%       'crop': [true | <false>] gain a few pixels by cropping aggresively
%       'scale': [true | <false>] export in my preferred size (w=768,h=576)
%
% NOTE: existing files will be overwritten!

%% Inputs
if nargin == 0, help ngraf.savehc, return, end
narginchk(1,inf)
p = inputParser();
p.addParameter('name','newfig',@(x)isstring(x)||(ischar(x)&&isrow(x)))
p.addParameter('crop',false,@(x)isscalar(x)&&islogical(x))
p.addParameter('scale',false,@(x)isscalar(x)&&islogical(x))
p.parse(varargin{:});
pr = p.Results;

%% matlab fig file
savefig(fh, pr.name, 'compact'); % compact requires >R2014a.

%% PDF and PNG
if pr.scale
    pos = fh.Position;
    fh.Position(3:4) = [768, 576];
    cu = onCleanup(@()set(fh,'position',pos));
end
if pr.crop
    export_fig(fh, pr.name, '-pdf', '-png', '-r600')
else
    export_fig(fh, pr.name, '-pdf', '-png', '-nocrop', '-r600')
end

%% BYU
end
