function [fh, ah, lh] = hist_nmoi(chain, varargin)
% Sample histogram of normalized moment of inertia values.

%% Input parsing
if nargin == 0
    fprintf('Usage:\n\thist_nmoi(chain, ''name'', ''value'')\n')
    return
end
narginchk(1,inf)
p = plots.basic_parser();
p.addParameter('norm','pdf')
p.addParameter('nbins',[],@(x)isscalar(x)&&x>0)
p.parse(varargin{:})
pr = p.Results;

%% Prepare the data
TOFS = [chain.tof];
ICE = double([TOFS.NMoI]);

%% Prepare the canvas
if isempty(pr.axes)
    [fh, ah] = ngraf.get_canvas(pr.target);
else
    ah = pr.axes;
    axes(ah)
    hold(ah, 'on')
end

%% Plot the histogram
if isempty(pr.nbins)
    lh = histogram(ICE, 'Normalization', pr.norm);
else
    lh = histogram(ICE, pr.nbins, 'Normalization', pr.norm);
end

%% Style and annotate axes
if isempty(pr.axes)
    xlabel('Normalized moment of inertia, $I/Ma_0^2$')
    ah.XLim = [min(ICE), max(ICE)];
    ah.XTick = ah.XLim;
    ah.YTickLabel = [];
    if ~isempty(pr.title)
        title(pr.title)
    end
end
end
