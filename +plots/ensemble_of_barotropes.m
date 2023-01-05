function lh = ensemble_of_barotropes(planets, varargin)
% Visualize ensemble barotropes with representative subset.

%% Input parsing
if nargin == 0
    fprintf('Usage:\n\tensemble_of_barotropes(planets, ''name'', ''value'')\n\t')
    fprintf('specify number of lines to plot with ''nlines'' (default 10)\n\t')
    fprintf('control transparency with ''alfa'' (default 1.0)\n')
    return
end
narginchk(1,inf)
p = plots.basic_parser();
p.addParameter('alfa',1.0,@(x)isscalar(x)&&x>0)
p.addParameter('nlines',10,@(x)isscalar(x)&&x>0)
p.parse(varargin{:})
pr = p.Results;
if pr.nlines > 500
    error('Too many lines can hang MATLAB (comment out this line to override)')
end
if isempty(pr.color), pr.color = 'byrc'; end

%% Prepare the data
pees = [planets.Pi]/1e11;
rhos = [planets.rhoi]/1000;
nsamples = size(pees,2);
rcs = rhos(end,:);
[~,I] = sort(rcs);
rhos = rhos(:,I);
pees = pees(:,I);

%% Prepare the canvas
if isempty(pr.axes)
    [~, ah] = ngraf.get_canvas(pr.target);
else
    ah = pr.axes;
    axes(ah)
end
hold(ah, 'on')
lh = gobjects(1,pr.nlines);

%% Plot the lines (pressure in Mbar density in 1000 kg/m^3)
map = pr.cmap(64);
skip = ceil(nsamples/pr.nlines);
L = 1:nsamples;
for k=1:skip:nsamples
    x = pees(:,k);
    y = rhos(:,k);
    lh(k) = line(x, y);
    
    % Apply color, significant or no
    if isequal(pr.color, 'byrc')
        d = (L(k) - min(L))/range(L); % hilariously inefficient
        cind = floor(d*length(map));
        cind = max(min(cind, 64), 1);
        lh(k).Color = map(cind,:);
    else
        lh(k).Color = pr.color;
    end
    try
        lh(k).Color(4) = pr.alfa;
    catch % 4-element color is undocumented; if fails just leave it
    end
end

%% Style and annotate axes
if isempty(pr.axes)
    ah.XLim(1) = 1e-5; % ignore very low pressure
    ah.YLim(1) = 1e-3;
    ah.XTick = 10.^(-5:1);
    ah.XScale = 'log';
    ah.YScale = 'log';
    xlabel('$P$ [Mbar]')
    ylabel('$\rho$ [1000 kg/m$^3$]')
    colormap(map)
    caxis([L(1),L(end)])
end
end
