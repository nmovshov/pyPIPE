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
p.addParameter('showad',false)
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
    hold(ah, 'on')
end

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
h = vline(min(pees(end,:)), '--');
h.Color = lh.FaceColor;

%% Reference adiabat overlays
if pr.showad
    hhe = barotropes.SCVH.HHE.isenP10T150Y275();
    h2o = barotropes.SCVH.ANEOSWATER.water_on_isenP10T150Y275();
    sio2 = barotropes.SCVH.ANEOSSERPENTINE.serpentine_on_isenP10T150Y275();
    iron = barotropes.SCVH.ANEOSIRON.iron_on_isenP10T150Y275();
    metalist = [10, 100];
    Z_mix = ahelpers.metal2z(metalist, 0.275);
    for k=1:length(metalist)
        mets(k) = barotropes.IdealMix([hhe,h2o],[1 - Z_mix(k),Z_mix(k)]); %#ok<AGROW>
    end
    p = logspace(6,13);
    rho_hhe = hhe.density(p);
    rho_h2o = h2o.density(p);
    for k=1:length(metalist)
        rho_met(k,:) = mets(k).density(p); %#ok<AGROW>
    end
    rho_sio2 = sio2.density(p);
    rho_iron = iron.density(p);

    y_hhe = rho_hhe/1000;
    y_h2o = rho_h2o/1000;
    y_met = rho_met/1000;
    y_sio2 = rho_sio2/1000;
    y_iron = rho_iron/1000;

    ah.YLimMode = 'manual';
    x = p/1e11;
    plot(x,y_hhe,'k-','DisplayName','H/He adiabat');
    lstyles = {'--',':'};
    for k=1:length(mets)
        plot(x,y_met(k,:),...
            Color='k',LineStyle=lstyles{k},LineWidth=2,...
            DisplayName=sprintf('H/He + H$_2$O, $%d\\times$ solar',metalist(k)));
    end
    plot(x,y_h2o,'k-.','DisplayName','H$_2$O');
    plot(x,y_sio2,'k-s','DisplayName','Serpentine');
    plot(x,y_iron,'k-d','DisplayName','Iron');
end

%% Misc
legend('Location','SE',FontSize=14)
