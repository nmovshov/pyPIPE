function [fh, ah, lh] = ensemble_of_barotropes(pees, rhos, varargin)
% Visualize ensemble barotropes with representative subset.

%% Input parsing
if nargin == 0
    fprintf('Usage:\n\tensemble_of_barotropes(pees,rhos, ''name'', ''value'')\n\t')
    fprintf('pees and rhos are nlevels-by-nsamples.\n\t')
    fprintf('control normalization with ''norm'' (default ''sqrt'')\n\t')
    fprintf('specify number of lines to plot with ''nlines'' (default 10)\n\t')
    fprintf('control transparency with ''alfa'' (default 1.0)\n')
    return
end
narginchk(1,inf)
p = plots.basic_parser();
p.addParameter('sqrtnorm',false)
p.addParameter('alfa',1.0,@(x)isscalar(x)&&x>0)
p.addParameter('nlines',10,@(x)isscalar(x)&&x>0)
p.parse(varargin{:})
pr = p.Results;
if pr.nlines > 500
    error('Too many lines can hang MATLAB (if you are sure comment out this line)')
end
nlevels = size(pees,1);
nsamples = size(pees,2);
if isempty(pr.color), pr.color = 'byrc'; end

%% Prepare the data
rcs = rhos(end,:);
[~,I] = sort(rcs);
rhos = rhos(:,I);
pees = pees(:,I);

%% Prepare the canvas
if isempty(pr.axes)
    [fh, ah] = ngraf.get_canvas(pr.target);
    lh = gobjects(1,pr.nlines);
else
    ah = pr.axes;
    axes(ah)
    hold(ah, 'on')
    fh = ah.Parent;
    lh = gobjects(1,pr.nlines);
end

%% Plot the lines (pressure in Mbar density in 1000 kg/m^3)
map = pr.cmap(64);
skip = ceil(nsamples/pr.nlines);
L = 1:nsamples;
for k=1:skip:nsamples
    x = pees(:,k)/1e11;
    y = rhos(:,k)/1000;
    if pr.sqrtnorm
        y = y./sqrt(x);
    end
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
    ah.XLim(1) = 1e-3; % ignore very low pressure
    ah.XTick = 10.^(-3:1);
    ah.XScale = 'log';
    ah.YScale = 'log';
    xlabel('$P$ [Mbar]')
    ylabel('$\rho$ [1000 kg/m$^3$]')
    if pr.sqrtnorm
        ah.YScale = 'linear';
        ylabel('$\rho \mathrm{[1000 kg/m^3]}/\sqrt{P}$')
    end
    colormap(map)
    caxis([L(1),L(end)])
    if ~isempty(pr.title)
        title(pr.title)
    end
end
end
