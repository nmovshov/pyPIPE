function [fh, ah, lh] = sigma_rho_of_r(profs, varargin)
% Level-by-level sample relative spread.

%% Input parsing
if nargin == 0
    fprintf('Usage:\n\tsigma_rho_of_r(profs,''name'',value)\n')
    fprintf('profs is nlevels-by-nsamples array\n')
    return
end
narginchk(1,inf)
p = plots.basic_parser();
p.addParameter('fractional',true,@(x)isscalar(x)&&islogical(x))
p.addParameter('svec',[],@(x)isvector(x))
p.parse(varargin{:})
pr = p.Results;
N = size(profs,1);
if isempty(pr.svec), pr.svec = linspace(1, 1/N, N)'; end

%% Prepare the data (assume chain holds fixed level-radii tofs!)
% These are the (assumed) common radii
s = pr.svec';
x = [(s/s(1)), 0];

% Get the median and +/- 1 and 2 sigma for emphasis
profs = profs';
y_med = median(profs);
y_med = [y_med, y_med(end)];
y_1sig = prctile(profs, [16, 84]);
y_1sig(:,end+1) = y_1sig(:,end);

% The quantity to plot is the width of the 1-sigma spread
y = y_1sig(2,:) - y_1sig(1,:);
if pr.fractional
    y = y./y_med; % fractional
else
    y = y/1000;   % or in 1000 kg/m^3
end

%% Prepare the canvas
if isempty(pr.axes)
    [fh, ah] = ngraf.get_canvas(pr.target);
else
    ah = pr.axes;
    fh = ah.Parent;
    axes(ah)
    hold(ah, 'on')
end

%% Plot the lines
lh = plot(x(2:end),y(2:end)); % first point is noise
lh.DisplayName = pr.label;

%% Style and annotate axes
if isempty(pr.axes)
    ah.YLim(1) = 0;
    xlabel('Level surface radius, $s/R_\mathrm{m}$')
    if pr.fractional
        ylabel('$\delta\rho/\widetilde{\rho}$')
    else
        ylabel('$\delta\rho$ [1000 kg/m$^3$]')
    end
    if ~isempty(pr.title)
        title(pr.title)
    end
end
end
