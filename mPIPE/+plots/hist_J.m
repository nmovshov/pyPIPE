function [fh, ah, lh] = hist_J(chain, n, varargin)
% Sample histogram of J_n values.

%% Input parsing
if nargin == 0
    fprintf('Usage:\n\thist_J(chain, n, ''name'', ''value'')\n')
    return
end
narginchk(2,inf)
p = plots.basic_parser();
p.addParameter('norm','pdf')
p.addParameter('nbins',[],@(x)isscalar(x)&&x>0)
p.parse(varargin{:})
pr = p.Results;

%% Prepare the data
TOFS = [chain.tof];
JAYS = double([TOFS.(sprintf('J%i',n))]);

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
    lh = histogram(JAYS, 'Normalization', pr.norm);
else
    lh = histogram(JAYS, pr.nbins, 'Normalization', pr.norm);
end

%% Style and annotate axes
if isempty(pr.axes)
    xlabel(sprintf('$J_%i$',n))
    ah.XLim = [min(JAYS), max(JAYS)];
    ah.XTick = ah.XLim;
    ah.YTickLabel = [];
    if ~isempty(pr.title)
        title(pr.title)
    end
end
end
