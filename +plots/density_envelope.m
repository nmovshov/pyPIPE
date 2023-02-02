function lh = density_envelope(planets, varargin)
% Level-by-level median and percentiles.

%% Input parsing
if nargin == 0
    fprintf('Usage:\n\taverage_of_profs(planets, ''name'', ''value'')\n\t')
    fprintf('set color with ''color''\n\t')
    fprintf('control transparency with ''alfa'' (default 1.0)\n')
    return
end
narginchk(1,inf)
p = plots.basic_parser();
p.addParameter('prctile',2);
p.addParameter('alfa',1.0,@(x)isscalar(x)&&x>0)
p.parse(varargin{:})
pr = p.Results;

%% Prepare the data (assume chain holds fixed level-radii tofs!)
x = planets(1).si/planets(1).si(1);
rhos = [planets.rhoi];
ylo = prctile(rhos, pr.prctile, 2)/1000;
yhi = prctile(rhos, 100 - pr.prctile, 2)/1000;

%% Prepare the canvas
if isempty(pr.axes)
    [~, ah] = ngraf.get_canvas(pr.target);
else
    ah = pr.axes;
    axes(ah)
end
hold(ah, 'on')

%% Plot the area (NOTE: density in 1000 kg/m^3)
[q,w,lh] = fill_between(x,ylo,yhi,[]);
if ~isempty(pr.color), lh.FaceColor = pr.color; end
lh.FaceAlpha = pr.alfa;
lh.LineStyle = 'none';
lh.DisplayName = pr.label;
delete([q,w]);

%% Style and annotate axes
if isempty(pr.axes)
    ah.YLim(1) = 0;
    xlabel('Level surface radius, $s/R_m$')
    ylabel('$\rho$ [1000 kg/m$^3$]')
    if ~isempty(pr.title)
        title(pr.title)
    end
end

%% Legend
if pr.legend
    legend(ah, 'off')
    legend(ah);
end
end
