function lh = barotrope_envelope(planets, varargin)
% Visualize sample envelope of rho(P).

%% Input parsing
if nargin == 0
    fprintf('Usage:\n\tbarotrope_envelope(planets, ''name'', ''value'')\n\t')
    fprintf('set color with ''color'' (default grey)\n\t')
    fprintf('control transparency with ''alfa'' (default 1.0)\n')
    return
end
narginchk(1,inf)
p = plots.basic_parser();
p.addParameter('adiabats',[])
p.addParameter('alfa',1.0,@(x)isscalar(x)&&x>0)
p.parse(varargin{:})
pr = p.Results;
interpolant = @spline;

%% Prepare the data
pees = [planets.Pi]/1e11;
rhos = [planets.rhoi]/1000;
x = logspace(-6, ceil(log10(max(pees(:)))), 1024)';
y = nan(length(x),size(rhos,2));
for k=1:size(rhos,2)
    ind = x < pees(end,k);
    y(ind,k) = interpolant(pees(:,k),rhos(:,k),x(ind));
end
[ylo, yhi] = ahelpers.density_prctiles(y, 2);
ylo(end) = []; yhi(end) = [];
ylo = ylo(:); yhi = yhi(:);

%% Prepare the canvas
if isempty(pr.axes)
    [~, ah] = ngraf.get_canvas(pr.target);
else
    ah = pr.axes;
    axes(ah);
end
hold(ah, 'on')

%% Plot the area (pressure in Mbar density in 1000 kg/m^3)
ind = ~isnan(ylo);
[q,w,lh] = fill_between(x(ind),ylo(ind),yhi(ind),[]);
lh.LineStyle = 'none';
if ~isempty(pr.color)
    lh.FaceColor = pr.color;
end
lh.FaceAlpha = pr.alfa;
lh.DisplayName = pr.label;
delete([q,w]);

%% Style and annotate axes
if isempty(pr.axes)
    ah.XLim(1) = 1e-5;
    ah.YLim(1) = 1e-3;
    ah.XTick = 10.^(-5:1);
    ah.XScale = 'log';
    ah.YScale = 'log';
    xlabel('$P$ [Mbar]')
    ylabel('$\rho$ [1000 kg/m$^3$]')
end
h = vline(min(pees(end,:)), '-');
h.Color = lh.FaceColor;
h.LineWidth = 0.5;

%% Reference adiabat overlays
for k=1:length(pr.adiabats)
    ad = pr.adiabats(k);
    p = logspace(6,13);
    rho = ad.density(p);
    y = rho/1000;
    ah.YLimMode = 'manual';
    x = p/1e11;
    lstyles = {'-','--',':','-.','-s','-d'};
    plot(x,y,color='k',LineStyle=lstyles{k},LineWidth=2,DisplayName=ad.name);
end

%% Misc
legend('Location','SE',FontSize=14)
