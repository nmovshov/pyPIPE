function [fh, ah, lh] = ensemble_of_profs(profs, varargin)
% Visualize ensemble with representative subset.

%% Input parsing
if nargin == 0
    fprintf('Usage:\n\tensemble_of_profs(profs, ''name'', ''value'')\n\t')
    fprintf('profs is nlevels-by-nsamples array\n\t')
    fprintf('specify number of lines to plot with ''nlines'' (default 100)\n\t')
    fprintf('control transparency with ''alfa'' (default 0.4)\n')
    return
end
narginchk(1,inf)
p = plots.basic_parser();
p.addParameter('alfa',0.4,@(x)isscalar(x)&&x>0)
p.addParameter('nlines',100,@(x)isscalar(x)&&x>0)
p.addParameter('svec',[],@(x)isvector(x))
p.parse(varargin{:})
pr = p.Results;
if pr.nlines > 500
    error('Too many lines can hang MATLAB (if you are sure comment out this line)')
end
N = size(profs,1);
if isempty(pr.color), pr.color = 'byrc'; end
if isempty(pr.svec), pr.svec = linspace(1, 1/N, N)'; end

%% Prepare the data
rcs = profs(end,:);
[~,I] = sort(rcs);
profs = profs(:,I);

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

%% Plot the lines (density in 1000 kg/m^3)
map = pr.cmap(64);
skip = ceil(size(profs,2)/pr.nlines);
L = 1:size(profs,2);
for k=1:skip:size(profs,2)
    x = [pr.svec/pr.svec(1); 0];
    y = [profs(:,k); profs(end,k)];
    lh(k) = line(x, y/1000);
    
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
    ah.YLim(1) = 0;
    xlabel('Level surface radius, $s/R_m$')
    ylabel('$\rho$ [1000 kg/m$^3$]')
    colormap(map)
    caxis([L(1),L(end)])
    if ~isempty(pr.title)
        title(pr.title)
    end
end
end
