function [fh, ah, lh] = polynomial_planets(filename, varargin)
% The ice giants as smooth polynomials.

%% Inputs
if nargin == 0
    fprintf('Usage:\n\tpolynomial_planets(filename,[''showgf''=true])\n')
    return
end
narginchk(1,inf)
p = plots.basic_parser();
p.addParameter('showgf',true)
p.parse(varargin{:})
pr = p.Results;

%% Load models
bmods = load(filename);
chain = bmods.pipe.chain;

%% Start with a clean canvas
[fh, ah] = ngraf.get_canvas(pr.target);
lh = gobjects(size(chain));

%% Make lines
for k=1:length(chain)
    link = chain(k);
    tof = link.tof;
    lh(k) = plot(tof.si/tof.s0, tof.rhoi/1000);
    lh(k).DisplayName = sprintf('degree %d',sum(isfinite(link.x))+1);
    lh(k).LineStyle = pr.lstyle;
    if ~isempty(pr.color), lh(k).Color = pr.color; end
    if pr.showgf
        lh(k).DisplayName = sprintf('%s; %0.2g-$\\sigma$ fit',...
            lh(k).DisplayName, mahal2sig(link.loss,3));
    end
end

%% Axis labels
ah.XLabel.String = 'Normalized mean radius $s/R_m$';
ah.YLabel.String = 'density [1000 kg/m$^3$]';

%% A legend
if pr.legend
    legend(lh)
end

%% Optional title
if isequal(pr.target, 'screen')
    title(bmods.pipe.name)
end

%% BYU
end
