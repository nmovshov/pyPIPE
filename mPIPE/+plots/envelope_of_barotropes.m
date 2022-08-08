function [fh, ah, h] = envelope_of_barotropes(pees, rhos, varargin)
% Visualize sample envelope of rho(P).

%% Input parsing
if nargin == 0
    fprintf('Usage:\n\tenvelope_of_barotropes(pees,rhos, ''name'', ''value'')\n\t')
    fprintf('pees and rhos are nlevels-by-nsamples.\n\t')
    fprintf('pressure in Mbar and density in 1000 kg/m^3.\n\t')
    fprintf('control normalization with ''sqrtnorm'' (default false)\n\t')
    fprintf('set color with ''color'' (default grey)\n\t')
    fprintf('control transparency with ''alfa'' (default 1.0)\n')
    return
end
narginchk(1,inf)
p = plots.basic_parser();
p.addParameter('sqrtnorm',false)
p.addParameter('alfa',1.0,@(x)isscalar(x)&&x>0)
p.addParameter('name','')
p.parse(varargin{:})
pr = p.Results;
interpolant = @spline;

%% Prepare the data
x = logspace(-6, ceil(log10(max(pees(:)))), 1024)';
y = nan(length(x),size(rhos,2));
for k=1:size(rhos,2)
    ind = x < pees(end,k);
    y(ind,k) = interpolant(pees(:,k),rhos(:,k),x(ind));
end
[ylo, yhi] = ahelpers.density_prctiles(y, 2);
ylo(end) = []; yhi(end) = [];
ylo = ylo(:); yhi = yhi(:);
if pr.sqrtnorm
    ylo = ylo./sqrt(x);
    yhi = yhi./sqrt(x);
end

%% Prepare the canvas
if isempty(pr.axes)
    [fh, ah] = ngraf.get_canvas(pr.target);
else
    ah = pr.axes;
    axes(ah)
    hold(ah, 'on')
    fh = ah.Parent;
end

%% Plot the area (pressure in Mbar density in 1000 kg/m^3)
ind = ~isnan(ylo);
[q,w,h] = fill_between(x(ind),ylo(ind),yhi(ind),[]);
h.LineStyle = 'none';
if ~isempty(pr.color)
    h.FaceColor = pr.color;
end
h.FaceAlpha = pr.alfa;
h.DisplayName = pr.name;
delete([q,w]);

%% Style and annotate axes
if isempty(pr.axes)
    ah.XLim(1) = 1e-3;
    ah.XTick = 10.^(-3:1);
    ah.XScale = 'log';
    ah.YScale = 'log';
    xlabel('$P$ [Mbar]')
    ylabel('$\rho$ [1000 kg/m$^3$]')
    if pr.sqrtnorm
        ah.YScale = 'linear';
        ylabel('$\rho\ \mathrm{[1000 kg/m^3]}/\sqrt{P}\ \mathrm{[Mbar]}$')
    end
end
