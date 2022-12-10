function [fh, ah, lh] = average_of_profs(profs, varargin)
% Level-by-level median and percentiles.

%% Input parsing
if nargin == 0
    fprintf('Usage:\n\taverage_of_profs(profs, ''name'', ''value'')\n\t')
    fprintf('profs is nlevels-by-nsamples array\n')
    return
end
narginchk(1,inf)
p = plots.basic_parser();
p.addParameter('shows1',true,@(x)isscalar(x)&&islogical(x))
p.addParameter('shows2',true,@(x)isscalar(x)&&islogical(x))
p.addParameter('s1color',[0.25,0.25,0.25],@(x)isrow(x)&&numel(x)==3)
p.addParameter('s2color',[0.50,0.50,0.50],@(x)isrow(x)&&numel(x)==3)
p.addParameter('s1alpha',0.75,@(x)isscalar(x)&&isreal(x))
p.addParameter('s2alpha',0.50,@(x)isscalar(x)&&isreal(x))
p.addParameter('svec',[],@(x)isvector(x))
p.parse(varargin{:})
pr = p.Results;
N = size(profs,1);
if isempty(pr.svec), pr.svec = linspace(1, 1/N, N)'; end

%% Prepare the data (assume chain holds fixed level-radii tofs!)
s = pr.svec';
rhos = profs';
x = [(s/s(1)), 0];
y = median(rhos);
y02 = prctile(rhos, 2);
y16 = prctile(rhos, 16);
y84 = prctile(rhos, 84);
y98 = prctile(rhos, 98);

y   = [y, y(end)];
y02 = [y02, y02(end)];
y16 = [y16, y16(end)];
y84 = [y84, y84(end)];
y98 = [y98, y98(end)];

%% Prepare the canvas
if isempty(pr.axes)
    [fh, ah] = ngraf.get_canvas(pr.target);
    lh = gobjects(1,3);
else
    ah = pr.axes;
    axes(ah)
    hold(ah, 'on')
end

%% Plot the lines (NOTE: density in 1000 kg/m^3)
lh(1) = plot(x,y/1000,'LineWidth',2);
lh(1).Color = 'k';
lh(1).DisplayName = 'median';

if pr.shows1
    [q,w,lh(2)] = fill_between(x,y16/1000,y84/1000,[],'facecolor',pr.s1color);
    lh(2).FaceAlpha = pr.s1alpha;
    lh(2).LineStyle = 'none';
    lh(2).DisplayName = '16th--84th prctile';
    delete([q,w]);
end
if pr.shows2
    [q,w,lh(3)] = fill_between(x,y02/1000,y98/1000,[],'facecolor',pr.s2color);
    lh(3).FaceAlpha = pr.s2alpha;
    lh(3).LineStyle = 'none';
    lh(3).DisplayName = '2nd--98th prctile';
    delete([q,w])
end
lh = lh(isgraphics(lh));

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
    legend(ah, flipud(lh), 'location','ne');
end
end
