function [fh, ah, lh] = average_of_profs_color(profs, varargin)
% Level-by-level rho(s) percentiles in color.

%% Input parsing
if nargin == 0
    fprintf('Usage:\n\taverage_of_profs_color(profs, ''name'', ''value'')\n\t')
    fprintf('profs is nlevels-by-nsamples array\n')
    return
end
narginchk(1,inf)
p = plots.basic_parser();
p.addParameter('clvls',10,@(x)isscalar(x)&&(x>0)&&(rem(x,1)==0))
p.addParameter('svec',[],@(x)isvector(x))
p.parse(varargin{:})
pr = p.Results;
N = size(profs,1);
if isempty(pr.svec), pr.svec = linspace(1, 1/N, N)'; end

%% Prepare the data
% These are the (assumed) common radii
s = pr.svec';
x = [(s/s(1)), 0];

% These are the percentiles we wish to colorize
prcs_lo = linspace(2,48,pr.clvls);
prcs_hi = 100 - prcs_lo;

% Get the data in the above percentiles
rhos = profs';
Y_LO = prctile(rhos, prcs_lo);
Y_HI = prctile(rhos, prcs_hi);
Y_LO(:,end+1) = Y_LO(:,end);
Y_HI(:,end+1) = Y_HI(:,end);

% Get the median and +/- 1 and 2 sigma for emphasis
y_med = median(rhos);
y_med = [y_med, y_med(end)];
y_1sig = prctile(rhos, [16, 84]);
y_2sig = prctile(rhos, [2, 98]);
y_1sig(:,end+1) = y_1sig(:,end);
y_2sig(:,end+1) = y_2sig(:,end);

%% Prepare the canvas
if isempty(pr.axes)
    [fh, ah] = ngraf.get_canvas(pr.target);
else
    ah = pr.axes;
    axes(ah)
    hold(ah, 'on')
end

%% Plot the lines (NOTE: density in 1000 kg/m^3)
lh = gobjects;

lh(end+1) = plot(x,y_med/1000,'LineWidth',2);
lh(end).Color = 'k';
lh(end).LineStyle = '-';
lh(end).DisplayName = 'median';
lh(end+1) = plot(x,y_1sig(1,:)/1000,'LineWidth',2);
lh(end).Color = 'k';
lh(end).LineStyle = '--';
lh(end).DisplayName = '1-$\sigma$';
lh(end+1) = plot(x,y_1sig(2,:)/1000,'LineWidth',2);
lh(end).Color = 'k';
lh(end).LineStyle = '--';
lh(end).DisplayName = '';
lh(end+1) = plot(x,y_2sig(1,:)/1000,'LineWidth',2);
lh(end).Color = 'k';
lh(end).LineStyle = ':';
lh(end).DisplayName = '2-$\sigma$';
lh(end+1) = plot(x,y_2sig(2,:)/1000,'LineWidth',2);
lh(end).Color = 'k';
lh(end).LineStyle = ':';
lh(end).DisplayName = '';

% Remove the placeholder that makes lh grow graphics instead of numbers
lh(1) = [];

% Plot the color zones with fill_between
cmap = pr.cmap(pr.clvls);
for k=1:pr.clvls
    [q,w,h] = fill_between(x,Y_LO(k,:)/1000,Y_HI(k,:)/1000,[],'facecolor',cmap(k,:));
    h.LineStyle = 'none';
    h.DisplayName = '';
    delete([q,w])
end
clear h

% Somehow the fill_between function gets the stacking of patches backwards
ah.Children(end-pr.clvls+1:end) = flipud(ah.Children(end-pr.clvls+1:end));

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
