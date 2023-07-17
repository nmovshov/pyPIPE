classdef samplooker
%% Some functions for looking at sampels of cooked planets
methods(Static)

    function hist_moi(planets, varargin)
        % Input
        p = samplooker.basic_parser();
        p.addParameter('rol',false,@(x)isscalar(x)&&islogical(x))
        p.parse(varargin{:})
        pr = p.Results;

        % Data
        ice = [planets.NMoI];
        mu = mean(ice); sig = std(ice);
        if pr.rol
            ice = ice(abs(ice - mu) < 3*sig);
        end

        % Canvas
        ah = samplooker.getx(pr.newfig, pr.target);

        % Plot
        lh = histogram(ice, pr.bins, 'Normalization', 'pdf');
        lh.DisplayName = pr.label;
        try
            validatestring('alfa',p.UsingDefaults);
        catch
            lh.FaceAlpha = pr.alfa;
        end
        try
            validatestring('color',p.UsingDefaults);
        catch
            lh.FaceColor = pr.color;
        end

        % Style and annotation
        xlabel('Normalized moment of inertia, $I/Ma_0^2$')
        ah.XTick = ah.XLim;
        ah.YTickLabel = [];
        if pr.legend
            legend(location='nw');
        end

    end

    function hist_rho0(planets, varargin)
        % Input
        p = samplooker.basic_parser();
        p.parse(varargin{:})
        pr = p.Results;

        % Data
        ros = [planets.rhoi];
        r0s = ros(1,:);

        % Canvas
        ah = samplooker.getx(pr.newfig, pr.target);

        % Plot
        lh = histogram(r0s, pr.bins, 'Normalization', 'pdf');
        lh.DisplayName = pr.label;
        try
            validatestring('alfa',p.UsingDefaults);
        catch
            lh.FaceAlpha = pr.alfa;
        end
        try
            validatestring('color',p.UsingDefaults);
        catch
            lh.FaceColor = pr.color;
        end

        % Style and annotation
        xlabel('1-bar density, $\rho_0$ [1000 kg/m$^3$]')
        ah.YTickLabel = [];
        if pr.legend
            legend(location='nw');
        end

    end

    function ensemble_of_profs(planets, varargin)
        % Input
        p = samplooker.basic_parser();
        p.addParameter('nlines',20)
        p.parse(varargin{:})
        pr = p.Results;

        % Data
        profs = [planets.rhoi];
        rcs = profs(end,:);
        [~,ind] = sort(rcs);
        profs = profs(:,ind);

        % Canvas
        ah = samplooker.getx(pr.newfig, pr.target); %#ok<NASGU>

        % Plot
        x = [planets(1).si/planets(1).si(1); 0];
        skip = ceil(length(planets)/pr.nlines);
        for k=1:skip:length(planets)
            y = [profs(:,k); profs(end,k)];
            lh = plot(x, y);
            try
                validatestring('color',p.UsingDefaults);
            catch
                lh.Color = pr.color;
            end
            lh.Color(4) = pr.alfa;
        end

        % Style and annotation
        xlabel('Level surface radius, $s/R_m$')
        ylabel('$\rho$ [1000 kg/m$^3$]')
    end

    function ensemble_of_barotropes(planets, varargin)
        % Input
        p = samplooker.basic_parser();
        p.addParameter('nlines',20)
        p.parse(varargin{:})
        pr = p.Results;

        % Data
        pees = [planets.Pi]*1e-11;
        rhos = [planets.rhoi]*1e-3;
        rcs = rhos(end,:);
        [~,ind] = sort(rcs);
        pees = pees(:,ind);
        rhos = rhos(:,ind);

        % Canvas
        ah = samplooker.getx(pr.newfig, pr.target);

        % Plot
        skip = ceil(length(planets)/pr.nlines);
        for k=1:skip:length(planets)
            x = pees(:,k);
            y = rhos(:,k);
            lh = plot(x, y);
            try
                validatestring('color',p.UsingDefaults);
            catch
                lh.Color = pr.color;
            end
            lh.Color(4) = pr.alfa;
        end

        % Style and annotation
        ah.XLim(1) = 1e-5;
        ah.YLim = [1e-3, 100];
        ah.XScale = 'log';
        ah.YScale = 'log';
        xlabel('$p$ [Mbar]')
        ylabel('$\rho$ [1000 kg/m$^3$]')
    end

    function density_envelope(planets, varargin)
        % Input
        p = samplooker.basic_parser();
        p.addParameter('prctile',2)
        p.parse(varargin{:})
        pr = p.Results;

        % Data
        profs = [planets.rhoi];
        rhos = profs';
        prcs_lo = pr.prctile;
        prcs_hi = 100 - prcs_lo;
        x = [planets(1).si/planets(1).si(1); 0];
        ylo = prctile(rhos, prcs_lo)/1000;
        yhi = prctile(rhos, prcs_hi)/1000;
        ylo = [ylo, ylo(end)]';
        yhi = [yhi, yhi(end)]';

        % Canvas
        ah = samplooker.getx(pr.newfig, pr.target); %#ok<NASGU>

        % Plot
        [q,w,lh] = fill_between(x,ylo,yhi,[]);
        lh.LineStyle = 'none';
        lh.DisplayName = pr.label;
        delete([q,w])
        try
            validatestring('color',p.UsingDefaults);
        catch
            lh.FaceColor = pr.color;
        end
        try
            validatestring('alfa',p.UsingDefaults);
        catch
            lh.FaceAlpha = pr.alfa;
        end

        % Style and annotation
        xlabel('Level surface radius, $s/R_m$')
        ylabel('$\rho$ [1000 kg/m$^3$]')
        if pr.legend && ~isempty(pr.label)
            legend(location='ne');
        end
    end

    function barotrope_envelope(planets, varargin)
        % Input
        p = samplooker.basic_parser();
        p.addParameter('prctile',2)
        p.parse(varargin{:})
        pr = p.Results;

        % Data
        pees = 1e-11*[planets.Pi];
        rhos = 1e-3*[planets.rhoi];
        x = logspace(-6, ceil(log10(max(max(pees)))), 1024)';
        y = nan(1024, length(planets));
        for k=1:length(planets)
            rhoofp = griddedInterpolant(pees(:,k), rhos(:,k),'pchip','none');
            y(:,k) = rhoofp(x);
        end
        prcs_lo = pr.prctile;
        prcs_hi = 100 - prcs_lo;
        ylo = prctile(y, prcs_lo, 2);
        yhi = prctile(y, prcs_hi, 2);

        % Canvas
        ah = samplooker.getx(pr.newfig, pr.target);

        % Plot
        ind = ~isnan(ylo);
        [q,w,lh] = fill_between(x(ind),ylo(ind),yhi(ind),[]);
        lh.LineStyle = 'none';
        lh.DisplayName = pr.label;
        delete([q,w])
        try
            validatestring('color',p.UsingDefaults);
        catch
            lh.FaceColor = pr.color;
        end
        try
            validatestring('alfa',p.UsingDefaults);
        catch
            lh.FaceAlpha = pr.alfa;
        end

        % Style and annotation
        ah.XScale = 'log';
        ah.YScale = 'log';
        ah.XLim(1) = 1e-5;
        ah.YLim = [1e-3, 100];
        xlabel('$p$ [Mbar]')
        ylabel('$\rho$ [1000 kg/m$^3$]')
        if pr.legend && ~isempty(pr.label)
            legend(location='ne');
        end
    end

    function p = basic_parser()
        % Many input parameters are common to all plots.
        p = inputParser();
        p.addParameter('newfig',false,@(x)isscalar(x)&&islogical(x))
        p.addParameter('bins',30)
        p.addParameter('target','proj',@(x)ischar(x)&&isrow(x))
        p.addParameter('color',[])
        p.addParameter('legend',true,@(x)isscalar(x)&&islogical(x))
        p.addParameter('label','',@(x)ischar(x)&&isrow(x))
        p.addParameter('alfa',0.6,@(x)isscalar(x)&&isreal(x))
    end

    function ah = getx(newfig,target)
        if newfig || isempty(get(groot,'CurrentFigure'))
            [~, ah] = ngraf.get_canvas(target);
        else
            ah = gca();
            axes(ah)
            hold(ah, 'on')
        end

    end

end % methods
end % class
