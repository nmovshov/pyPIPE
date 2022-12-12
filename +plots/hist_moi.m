function lh = hist_moi(planets, varargin)
% Sample histogram of normalized moment of inertia values.

%% Input parsing
if nargin == 0
    fprintf('Usage:\n\thist_nmoi(planets, key=value,...)\n')
    return
end
narginchk(1,inf)
p = plots.basic_parser();
p.addParameter('showmu',false,@(x)isscalar(x)&&islogical(x))
p.addParameter('rol',false,@(x)isscalar(x)&&islogical(x))
p.addParameter('norm','pdf')
p.addParameter('nbins',[],@(x)isscalar(x)&&x>0)
p.parse(varargin{:})
pr = p.Results;

%% Prepare the data
ice = [planets.NMoI];
mu = mean(ice); sig = std(ice);
if pr.rol
    ice = ice(abs(ice - mu) < 3*sig);
end

%% Prepare the canvas
if isempty(pr.axes)
    [~, ah] = ngraf.get_canvas(pr.target);
else
    ah = pr.axes;
    axes(ah)
    hold(ah, 'on')
end

%% Plot the histogram
if isempty(pr.nbins)
    lh = histogram(ice, 'Normalization', pr.norm);
else
    lh = histogram(ice, pr.nbins, 'Normalization', pr.norm);
end
if ~isempty(pr.label)
    lh.DisplayName = pr.label;
end
if pr.showmu
    vline(mu,'r--');
end

%% Style and annotate
if isempty(pr.axes)
    xlabel('Normalized moment of inertia, $I/Ma_0^2$')
    ah.XLim = [min(ice), max(ice)];
    ah.XTick = ah.XLim;
    if pr.showmu
        ah.XTick = sort([ah.XTick, mu]);
    end
    ah.YTickLabel = [];
    if ~isempty(pr.title)
        title(pr.title)
    end
    if pr.legend
        legend(location='nw');
    end
end
end
