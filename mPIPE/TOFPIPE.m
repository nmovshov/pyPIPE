classdef TOFPIPE < handle
    %TOFPIPE TOF-based planetary-interior posterior explorer.
    
    %% Properties
    properties
        name
        chain
        opts
    end
    properties % optional
        pmodel
        lossfunction
        likefunction
        observables
        toflevels
        tofopts
    end
    properties (Dependent)
        nlinks
        multlinks
        badlinks
    end
    
    %% Simple constructor
    methods
        function obj = TOFPIPE(varargin)
            % The constructor only populates the options struct.
            obj.opts = pipeset(varargin{:});
        end
    end
    
    %% Public methods
    methods
        function add_link(obj, tof)
            new_link.x            = [];
            new_link.tof          = [];
            new_link.runtime      = [];
            new_link.loss         = [];
            new_link.like         = [];
            new_link.multiplicity = 1;
            
            if isempty(obj.chain)
                obj.chain = new_link;
            else
                obj.chain(end+1) = new_link;
            end
            
            if nargin > 1
                assert(isa(tof,'TOFPlanet') || isstruct(tof),...
                    'expected a TOFPlanet or reduced TOFPlanet here.')
                obj.chain(end).tof = tof;
                obj.chain(end).runtime = 0;
            end
        end
        
        function run_parametric_sequence(obj, X)
            % Run the sequence of parameters in array X through the PIPE.
            %
            % Usage: pipe.run_parameteric_sequence(X)
            %    X - nObservations-by-ndims array of parameters
            
            if nargin < 2
                help TOFPIPE.run_parametric_sequence
                return
            end
            if isempty(obj.pmodel)
                error('Set PIPE.pmodel field.')
            end
            
            for k=1:size(X,1)
                x = X(k,:);
                fprintf('Running model %g of %g...',k,size(X,1))
                t_seq = tic;
                if (k > 1) && isequal(x, X(k-1,:))
                    obj.chain(end).multiplicity = obj.chain(end).multiplicity + 1;
                else
                    obj.run_parametric_model(obj.pmodel, x);
                end
                t_seq = toc(t_seq);
                fprintf('done. (%g sec.)\n',t_seq)
            end
        end
        
        function run_parametric_model(obj, mdl, x)
            % Construct, run, and save a realization of a parametric model.
            %
            % Usage: pipe.run_parametric_model(mdl, x)
            %    mdl : function handle
            %      mdl(x) should return a TOFPlanet object ready to relax to HE.
            %    x : vector
            %      Inputs to mdl.
            
            % These are the inputs
            if nargin < 3
                help CMSPIPE.run_parametric_model
                return
            end
            narginchk(3,3)
            validateattributes(mdl,{'function_handle'},{},'','mdl',1)
            validateattributes(x,{'numeric'},{'row'},2)
            
            % Make sure we have at least default tofopts
            if isempty(obj.tofopts)
                obj.tofopts = tofset;
            end
            
            % And save the name of the mdl handle, just for fun
            if isempty(obj.pmodel)
                obj.pmodel = char(mdl);
            end
            
            % Optional communication
            verb = obj.opts.verbosity;
            if verb > 3
                fprintf('\n  Running %s(',char(mdl))
                for k=1:length(x), fprintf('%g,',x(k)); end
                fprintf('\b)...')
            end
            
            % Extend the chain to hold the new model
            obj.add_link;
            obj.chain(end).x = x;
            
            % Construct the model from given parameters and store its handle
            tof = realize_model(mdl, x, obj.tofopts);
            obj.chain(end).tof = tof;
            
            % Try to use last link to speed up convergence to HE
            if ~isempty(tof)
                if obj.opts.linkJ0 && obj.nlinks > 1
                    try
                        tof.set_ss_guesses(obj.chain(end-1).tof.ss);
                    catch ME
                        warning(ME.identifier,...
                            'Linking Ss failed with message:\n%s',ME.message)
                    end
                end
            end
            
            % Attempt to quick-reject model without converging
            if ~isempty(tof) && obj.opts.quickreject
                qr = quick_reject_model(tof);
            else
                qr = 0;
            end
            
            % Run appropriate TOF method to converge the model
            if ~isempty(tof) && (qr == 0)
                obj.chain(end).runtime = converge_model(tof, obj.opts);
            else
                obj.chain(end).runtime = NaN;
            end
            
            % Store loss or likelihood value if available
            if ~isempty(tof) && ~isempty(obj.lossfunction)
                obj.chain(end).loss = obj.lossfunction(tof);
            end
            if ~isempty(tof) && ~isempty(obj.likefunction)
                obj.chain(end).like = obj.likefunction(tof);
            end
            
            % Optional communication
            if verb > 3
                fprintf('done. (%g sec)\n', max(0, obj.chain(end).runtime))
                if isnan(obj.chain(end).runtime)
                    fprintf('  Model failed or was quick-rejected.\n')
                else
                    fprintf('  Putting converged model in link #%d.\n',obj.nlinks)
                end
            end
            
            % Optionally reduce link to save memory
            if ~isempty(tof) && (obj.opts.reducelinks)
                s = tof.to_struct(obj.opts.reducelinks, obj.opts.linkJ0);
                obj.chain(end).tof = s;
                tof.delete();
                if obj.nlinks > 1
                    obj.chain(end-1).tof.ss = [];
                end
            end
            
            % Optionally dump backup file as chain grows
            if ~mod(obj.nlinks, obj.opts.backupfreq)
                obj.dump_chain;
            end
            
            % ABYU
            return
        end
        
        function L = calculate_loss(obj, mdl, x, lfh, pen)
            % Run parametric model and return value of a loss function.
            %
            % Pass this function to optimization routines, e.g., fminsearch. The
            % loss function value will include a finite penalty for non-physical
            % parameters to help guide the optimization.
            %
            % Usage: pipe.calculate_loss(mdl, x, lfh)
            %    lfh : function handle
            %      loss function; lfh(tof) returns a positive number (some sort of
            %      distance from nominal model, see +losses package help)
            %    mdl, x : inputs to TOFPIPE.run_parametric_model
            
            if nargin < 4
                help TOFPIPE.calculate_loss
                return
            end
            if nargin < 5, pen = true; end % undocumented disable penalty to debug
            
            % Start by letting run_parametric do its thing
            obj.run_parametric_model(mdl, x);
            tof = obj.chain(end).tof;
            
            % Call the loss function, if there is a converged model
            if isempty(tof), L = inf; return, end
            L = lfh(tof);
            
            % Add severe but finite penalty for non-physical density
            if pen
                if any(tof.rhoi < 0)
                    yp = abs(sum(tof.rhoi(tof.rhoi < 0))/tof.rhobar);
                    yp = yp/numel(tof.rhoi);
                else
                    deltas = [tof.rhoi(1); diff(tof.rhoi(:))]/tof.rhobar;
                    ind = deltas < 0;
                    yp = abs(sum(deltas(ind))); % remember sum([]) is zero
                    yp = yp/numel(tof.rhoi);
                end
                yp = yp*1e6;
            else
                yp = 0;
            end
            
            % ABYU
            L = L + yp;
        end
        
        function [link, ind] = winner(obj)
            % Return link with lowest loss function value.
            
            L = [obj.chain.loss];
            [~, ind] = min(L);
            link = obj.chain(ind);
        end
        
        function [link, ind] = loser(obj)
            % Return link with highest loss function value.
            
            L = [obj.chain.loss];
            [~, ind] = max(L);
            link = obj.chain(ind);
        end
        
        function L = lossify_chain(obj, lfh)
            % Run a loss function evaluation on all chain links and return values.
            %
            % Usage: L = pipe.lossify_chain(lfh)
            %    lfh : function handle
            %      lfh(tof) returns loss/distance value (tip: help losses)
            
            if nargin < 2
                help TOFPIPE.lossify_chain
                return
            end
            narginchk(2,2)
            validateattributes(lfh,{'function_handle'},{},'','lfh',1)
            
            % Run on a clean chain
            obj.remove_bad_links;
            if isempty(obj.chain), L = []; return, end
            
            L = NaN(1,obj.nlinks);
            textprogressbar('Lossifying chain: ')
            cob = onCleanup(@()clear('textprogressbar'));
            for k=1:length(L)
                L(k) = lfh(obj.chain(k).tof);
                obj.chain(k).loss = L(k);
                obj.chain(k).like = -0.5*L(k)^2;
                textprogressbar(100*k/obj.nlinks)
            end
            textprogressbar(' done.')
            obj.lossfunction = lfh;
        end
        
        function ic = remove_bad_links(obj)
            % Remove links that failed to converge to physical a model.
            
            if isempty(obj.chain), return, end
            ind = cellfun(@(x)isempty(x)||isnan(x), {obj.chain.runtime});
            ic = sum(ind);
            obj.chain(ind) = [];
        end
        
        function ic = remove_losers(obj, thresh)
            % Remove links with loss function value above threshold.
            %
            % Usage: pipe.remove_losers(threshold)
            
            if nargin < 2
                help TOFPIPE.remove_losers
                return
            end
            validateattributes(thresh,{'numeric'},{'positive','scalar'})
            
            L = [obj.chain.loss];
            ind = (L > thresh);
            obj.chain(ind) = [];
            ic = sum(ind);
        end
        
        function reduce_chain(obj,lvl)
            % Reduce memory and disk usage by converting to static single structs.
            
            if nargin > 1, obj.opts.reducelinks = lvl; end
            obj.remove_bad_links;
            textprogressbar('Reducing chain: ')
            cob = onCleanup(@()clear('textprogressbar'));
            ochain = obj.chain;
            for k=1:length(ochain)
                if isa(ochain(k).tof, 'TOFPlanet')
                    s = ochain(k).tof.to_struct(obj.opts.reducelinks,false);
                    ochain(k).tof.delete(); % maybe helps garbage collection?
                else
                    s = structfun(@single, ochain(k).tof, 'UniformOutput', false);
                    s.ss = [];
                    s.name = obj.name;
                    if obj.opts.reducelinks > 2
                        s.si = [];
                        s.ai = [];
                        s.rhoi = [];
                        s.Pi = [];
                        s.mi = [];
                        s.zi = [];
                    end
                end
                ochain(k).tof = s;
                textprogressbar(100*k/obj.nlinks)
            end
            obj.chain = ochain;
            textprogressbar(' done.')
        end
        
        function echain = export_chain(obj)
            % Export cleaned up, flattened, and simplified chain.
            
            obj.remove_bad_links;
            textprogressbar('Exporting chain: ')
            cob = onCleanup(@()clear('textprogressbar'));
            n = obj.nlinks;
            ochain = obj.chain;
            echain = ochain(1);
            for k=2:n
                for j=1:ochain(k).multiplicity
                    echain(end+1) = ochain(k); %#ok<AGROW>
                end
                textprogressbar(100*k/n)
            end
            textprogressbar(' done.')
            if ~isempty(echain)
                echain = rmfield(echain, {'runtime', 'multiplicity'});
            end
        end
        
        function flush_chain(obj,pbar)
            % Clean up and flatten chain.
            
            if nargin < 2, pbar = true; end
            obj.remove_bad_links;
            if pbar
                textprogressbar('Flushing chain: ')
                cob = onCleanup(@()clear('textprogressbar'));
            end
            n = obj.nlinks;
            ochain = obj.chain;
            echain = ochain(1);
            for k=2:n
                for j=1:ochain(k).multiplicity
                    echain(end+1) = ochain(k); %#ok<AGROW>
                    echain(end).multiplicity = 1;
                end
                if pbar
                    textprogressbar(100*k/n)
                end
            end
            if pbar
                textprogressbar(' done.')
            end
            obj.chain = echain;
        end
        
        function V = getmdlfield(obj, tofprop)
            % Extract a field from chain's TOFs and return as vector or matrix.
            %
            % Usage: V = pipe.getmdlfield(tofprop)
            %    tofprop - a property of the saved tof structs.
            %
            % Example: M = pipe.getmdlfield('M')
            
            if nargin == 1
                help TOFPIPE.getmdlfield
                V = [];
                return
            end
            if isempty(obj.chain)
                warning('Empty chain.')
                V = [];
                return
            end
            narginchk(2,2)
            validateattributes(tofprop, {'char'}, {'row'})
            TOFS = [obj.chain.tof];
            try
                V = [TOFS.(tofprop)];
            catch ME
                help TOFPIPE.getmdlfield
                rethrow(ME)
            end
        end
        
        function X = getx(obj)
            % Retrieve obj.chain.x in single array.
            
            if obj.nlinks == 0
                X = [];
            else
                X = nan(length(obj.chain(1).x), obj.nlinks);
                for k=1:obj.nlinks
                    X(:,k) = obj.chain(k).x;
                end
            end
        end
        
    end % End of public methods block
    
    %% Visualizers
    methods
        function [fh, ah, lh] = plot_chain(obj, varargin)
            % Quick glance at ALL profiles.
            
            % Input parsing
            if obj.nlinks == 0, return, end
            if obj.nlinks > 30
                error('This is a long chain, did you mean to call plot_ensemble?')
            end
            p = basic_parser();
            p.addParameter('alfa',false,@(x)islogical(x)&&isscalar(x))
            p.parse(varargin{:})
            pr = p.Results;
            
            % Prepare the data
            links = obj.chain;
            L = 1:obj.nlinks;
            
            % Prepare the canvas
            if isempty(pr.axes)
                [fh, ah] = get_canvas(pr.target);
                lh = gobjects(size(links));
            else
                ah = pr.axes;
                axes(ah)
                hold(ah, 'on')
            end
            
            % Plot the lines (density in 1000 kg/m^3)
            map = pr.cmap(64);
            for k=1:length(L)
                tof = links(k).tof;
                x = [tof.si/tof.s0; 0];
                y = [tof.rhoi; tof.rhoi(end)];
                lh(k) = line(x, y/1000);
                
                % Apply color, significant or no
                d = (L(k) - min(L))/range(L); % hilariously inefficient
                cind = floor(d*length(map));
                cind = max(min(cind, 64), 1);
                lh(k).Color = map(cind,:);
                
                % Help out by partly alphaing successive lines (undocumented)
                if pr.alfa
                    lh(k).Color(4) = 1 - 0.8*k/length(L);
                end
            end
            
            % Style and annotate axes
            if isempty(pr.axes)
                ah.YLim(1) = 0;
                xlabel('Level surface radius, $s/R_m$')
                ylabel('$\rho$ [1000 kg/m$^3$]')
                colormap(map)
                if ~isscalar(L), caxis([L(1),L(end)]), end
            end
        end
        
        function [fh, ah, lh] = plot_ensemble_rho_of_r(obj, varargin)
            % Visualize ensemble with representative subset.
            
            % Input parsing
            if obj.nlinks == 0, return, end
            p = basic_parser();
            p.addParameter('alfa',0.4,@(x)isscalar(x)&&x>0)
            p.addParameter('nlines',100,@(x)isscalar(x)&&x>0)
            p.parse(varargin{:})
            pr = p.Results;
            if pr.nlines > 500
                error('Too many lines can hang MATLAB (if you are sure comment out this line)')
            end
            if isempty(pr.color), pr.color = 'byrc'; end
            
            % Prepare the data
            links = obj.chain;
            TOFS = [links.tof];
            rhos = [TOFS.rhoi];
            rcs = rhos(end,:);
            [~,I] = sort(rcs);
            links = links(I);
            
            % Prepare the canvas
            if isempty(pr.axes)
                [fh, ah] = get_canvas(pr.target);
                lh = gobjects(1,pr.nlines);
            else
                ah = pr.axes;
                axes(ah)
                hold(ah, 'on')
            end
            
            % Plot the lines (density in 1000 kg/m^3)
            map = pr.cmap(64);
            skip = ceil(length(links)/pr.nlines);
            L = 1:length(links);
            for k=1:skip:length(links)
                tof = links(k).tof;
                x = [tof.si/tof.s0; 0];
                y = [tof.rhoi; tof.rhoi(end)];
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
            
            % Style and annotate axes
            if isempty(pr.axes)
                ah.YLim(1) = 0;
                xlabel('Level surface radius, $s/R_m$')
                ylabel('$\rho$ [1000 kg/m$^3$]')
                colormap(map)
                if ~isscalar(L), caxis([L(1),L(end)]), end
            end
        end
        
        function [fh, ah, lh] = plot_averaged_rho_of_r(obj, varargin)
            % Level-by-level median and percentiles.
            
            % Input parsing
            if obj.nlinks == 0, return, end
            if obj.nlinks == 1
                warning('Single-link chain; call plot_chain() instead.')
                return
            end
            if any([obj.chain.multiplicity] > 1)
                error('Multiplicity detected; call OBJ.flush_chain to resolve multiple links.')
            end
            p = basic_parser();
            p.addParameter('shows1',true,@(x)isscalar(x)&&islogical(x))
            p.addParameter('shows2',true,@(x)isscalar(x)&&islogical(x))
            p.addParameter('s1color',[0.25,0.25,0.25],@(x)isrow(x)&&numel(x)==3)
            p.addParameter('s2color',[0.50,0.50,0.50],@(x)isrow(x)&&numel(x)==3)
            p.addParameter('s1alpha',0.75,@(x)isscalar(x)&&isreal(x))
            p.addParameter('s2alpha',0.50,@(x)isscalar(x)&&isreal(x))
            p.parse(varargin{:})
            pr = p.Results;
            
            % Prepare the data (assume chain holds fixed level-radii tofs!)
            s = obj.chain(1).tof.si';
            TOFS = [obj.chain.tof];
            rhos = [TOFS.rhoi]';
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
            
            % Prepare the canvas
            if isempty(pr.axes)
                [fh, ah] = get_canvas(pr.target);
                lh = gobjects(1,3);
            else
                ah = pr.axes;
                axes(ah)
                hold(ah, 'on')
            end
            
            % Plot the lines (NOTE: density in 1000 kg/m^3)
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
            
            % Style and annotate axes
            if isempty(pr.axes)
                ah.YLim(1) = 0;
                xlabel('Level surface radius, $s/R_m$')
                ylabel('$\rho$ [1000 kg/m$^3$]')
            end
            
            % Legend
            legend(ah, 'off')
            legend(ah, flipud(lh), 'location','ne');
        end
    end % End visualizers
    
    %% Private methods
    methods (Access = private)
        function dump_chain(pipe)
            persistent d
            if isempty(d), d = 0; end
            fname = sprintf('dump%d.mat',mod(d,2));
            try
                fprintf('Saving pipe to:  %s\n',fullfile(pwd,fname))
                save(fname,'pipe')
            catch ME
                warning(ME.identifier,'Backup dump failed with message:\n%s',ME.message)
            end
            d = d + 1;
        end
    end % End of private methods block
    
    %% Graveyard (obsolete methods)
    methods (Access = private)
        function accept = mhsample(obj, x0, nsamples, fun, prop)
            % Draw Metropolis-Hastings sample using an objective function.
            %
            % Usage:
            %   accept = pipe.mhsample(x0, nsamples, fun, prop)
            %     x0 - array of seed parameters
            %     nsamples - number of samples to draw
            %     fun - handle to objective function (returns LOSS value)
            %     prop - handle to proposal function y=f(x) (must be symmetric)
            %     accept - sample acceptance rate
            
            % These are the inputs
            if nargin == 1 && nargout == 0
                help TOFPIPE.mhsample
                return
            end
            narginchk(5,5)
            validateattributes(x0,{'numeric'},{'2d'},'','x0',1)
            validateattributes(nsamples,{'numeric'},{'integer','positive','scalar'},'','nsamples',2)
            validateattributes(fun,{'function_handle'},{},'','fun',3)
            validateattributes(prop,{'function_handle'},{},'','prop',4)
            
            % Local variables
            logify = @(y)-0.5*y^2; % loss-to-loglike
            verb = obj.opts.verbosity;
            
            % Begin by running the seed model
            if verb > 0
                fprintf('\nEvaluating seed point... \n')
            end
            t = tic;
            [L0, out0] = fun(x0);
            t = toc(t);
            L0 = logify(L0);
            if (L0 == -Inf)
                error('Sampling must start from finite likelihood seed point.')
            end
            if verb > 0
                fprintf('Evaluating seed point... done. (%g sec)\n', t)
                fprintf('Seed point scored L0 = %g.\n', L0)
            end
            obj.add_link();
            obj.chain(end) = out0.link;
            obj.chain(end).like = L0;
            if (obj.opts.reducelinks) % Optionally reduce tof to struct
                tof = out0.tof;
                s = tof.to_struct(obj.opts.reducelinks,false);
                obj.chain(end).tof = s;
                tof.delete();
            end
            
            % Now the Metropolis-Hastings algorithm
            t_smpl = tic;
            x = x0;
            Lx = L0;
            accept = 0;
            for n=1:nsamples
                y = prop(x);
                if verb > 0
                    fprintf('\nEvaluating proposal %d (of %d)... \n', n, nsamples)
                end
                t = tic;
                [Ly, outy] = fun(y);
                t = toc(t);
                Ly = logify(Ly);
                if verb > 0
                    fprintf('Evaluating proposal %d (of %d)... done. (%g sec)\n',n,nsamples,t)
                    fprintf('Proposal evaluation: Ly = %g (Lx = %g).\n', Ly, Lx)
                end
                acc = log(rand) < min(Ly - Lx, 0);
                if acc
                    x = y;
                    Lx = Ly;
                    obj.add_link();
                    obj.chain(end) = outy.link;
                    obj.chain(end).like = Ly;
                    if verb > 0
                        fprintf('Proposal accepted. ')
                    end
                else
                    obj.chain(end).multiplicity = obj.chain(end).multiplicity + 1;
                    if verb > 0
                        fprintf('Proposal rejected. ')
                    end
                end
                
                accept = accept+(acc);
                if verb > 0
                    fprintf('Running acceptance rate is %g.\n',accept/n)
                end
                
                % Optionally reduce tof to struct to save memory
                if (obj.opts.reducelinks) && acc
                    tof = outy.tof;
                    s = tof.to_struct(obj.opts.reducelinks,false);
                    obj.chain(end).tof = s;
                    tof.delete();
                end
                
                % Optionally dump backup file as sample grows
                if ~mod(n, obj.opts.backupfreq)
                    obj.dump_chain;
                end
            end
            
            % Return acceptance rate
            accept = accept/nsamples;
            ET = toc(t_smpl);
            if verb > 0
                fprintf('\nAll samples drawn. Acceptance rate was %g.\n',accept);
                fprintf('Run time %s\n',lower(seconds2human(ET)))
            end
        end
        
        function hitrate = essample(obj, nsamples, meantof, seedtof, SIG)
            % Draw layer densities by elliptical slice sampling around mean model.
            %
            % Usage:
            %   hitrate = pipe.essample(nsamples, meantof, seedtof, SIG)
            %     nsamples - number of samples to draw
            %     meantof - prior mean model, a converged TOFPlanet (or struct)
            %     seedtof - seed model, a converged TOFPlanet (or struct)
            %     SIG - prior covariance matrix on densities (help gpkernels)
            %   output:
            %     hitrate - average slice hit rate (1/bracket-shrinkings)
            
            % These are the inputs
            try
                narginchk(5,5)
                validateattributes(nsamples,...
                    {'numeric'},{'integer','nonnegative','scalar'},'','nsamples',1)
                validateattributes(meantof,...
                    {'TOFPlanet','struct'},{'scalar'},'','meantof',2)
                validateattributes(seedtof,...
                    {'TOFPlanet','struct'},{'scalar'},'','seedtof',3)
                N1 = length(meantof.si); N2 = length(seedtof.si);
                validateattributes(SIG,...
                    {'numeric'},{'2d','size',[N1,N2]},'','SIG',4)
                assert(norm(meantof.si/meantof.s0-seedtof.si/seedtof.s0)/N1 < 1e-3,...
                    'Mean model and seed model must share radius-grid (tof.si)')
                if isempty(obj.observables)
                    error('Set pipe.observables field (tip: help observables)')
                end
                if isempty(obj.likefunction)
                    error('Set pipe.likefunction field (tip: help loglikes)')
                end
            catch ME
                if (nargin == 1) && (nargout == 0)
                    help TOFPIPE.essample
                    return
                else
                    rethrow(ME)
                end
            end
            
            % We will use these local names a lot
            myobs = obj.observables;
            myN   = length(meantof.si);
            mymdl = @(~,x)tofmodels.rhoofs(x);
            mylogpdf = obj.likefunction;
            verb = obj.opts.verbosity;
            tof0 = meantof;
            rho0 = tof0.rhobar;
            y0 = [tof0.rhoi(1); diff(tof0.rhoi)]/rho0;
            assert(all(y0 >= 0), 'Mean model cannot have density inversions')
            f0 = log(y0);
            x0 = tof0.si;
            nbackup = obj.opts.backupfreq;
            obj.opts.backupfreq = 0; % override meaning of backupfreq
            
            % Allow dry run with nsamples==0
            if nsamples == 0
                fprintf('Dry run successful.\n')
                hitrate = 0;
                return
            end
            
            % Begin by running the mean model
            if verb > 0
                fprintf('\nEvaluating mean model...')
            end
            obj.run_parametric_model(mymdl, myN, [x0, tof0.rhoi], myobs, obj.tofopts);
            if isnan(obj.chain(end).runtime)
                error('Sampling requires a valid mean model.')
            end
            L0 = mylogpdf(obj.chain(end).tof, myobs);
            if ~isfinite(L0)
                error('Sampling requires a finite likelihood mean model.')
            end
            if verb > 0
                fprintf('done. (%g sec)\n', max(0, obj.chain(end).runtime))
                fprintf('Mean model scored L0 = %g.\n', L0)
            end
            
            % Then the seed model
            if verb > 0
                fprintf('\nEvaluating seed model...')
            end
            obj.run_parametric_model(mymdl, myN, [x0, seedtof.rhoi], myobs, obj.tofopts);
            if isnan(obj.chain(end).runtime)
                error('Sampling must start with a valid seed model.')
            end
            L0 = mylogpdf(obj.chain(end).tof, myobs);
            if ~isfinite(L0)
                error('Sampling must start from finite likelihood seed model.')
            end
            if verb > 0
                fprintf('done. (%g sec)\n', max(0, obj.chain(end).runtime))
                fprintf('Seed model scored L0 = %g.\n', L0)
            end
            
            % Now, the elliptical slice algorithm
            warning off TOFPIPE:realize_model_fail
            t_smpl = tic;
            scount = ones(1,nsamples);
            for n=1:nsamples
                if verb > 0
                    fprintf('\nDrawing sample %d (of %d)...\n', n, nsamples)
                end
                
                % Input: current state, f, from rhoi of last link.
                tof = obj.chain(end).tof;
                y = [tof.rhoi(1); diff(tof.rhoi)]/rho0;
                f = log(y) - f0;
                f(~isfinite(f)) = 0; % where inf arithmetics fails
                
                % Ellipse is defined by three points: 0, f, and nu
                nu = mvnrnd(zeros(myN,1),SIG)';
                nu(1) = 0; % respect boundary condition (surface density)
                
                % Set a log-likelihood threshold
                Ly = mylogpdf(tof, myobs) + log(rand);
                
                % Draw initial proposal from the full ellipse
                teta = rand*2*pi;
                [tetamin, tetamax] = deal(teta-2*pi, teta);
                
                % Shrink bracket until proposal is drawn on the L > Ly "slice"
                success = false;
                while (tetamax - tetamin) > eps
                    if verb > 1
                        fprintf('  Evaluating proposal (theta = %4.2f pi)...',...
                            teta/pi)
                    end
                    fp = f*cos(teta) + nu*sin(teta);
                    yp = exp(fp + f0);
                    yp(1) = y0(1); % respect boundary condition
                    rop = cumsum(yp)*rho0;
                    obj.run_parametric_model(mymdl, myN, [x0,rop], myobs, obj.tofopts);
                    tof = obj.chain(end).tof;
                    if isempty(tof) || isnan(obj.chain(end).runtime)
                        Lfp = -Inf;
                    else
                        Lfp = mylogpdf(tof, myobs);
                    end
                    if verb > 1
                        fprintf('done. Proposal scored L = %g (Ly = %g).\n',...
                            Lfp, Ly)
                    end
                    if Lfp > Ly
                        if verb > 1
                            fprintf('  Proposal on slice.\n')
                        end
                        success = true;
                        obj.chain(end).x = []; % waste of space...
                        obj.chain(end).x = teta; % interesting diagnostic
                        break
                    else
                        if verb > 1
                            fprintf('  Proposal off slice, shrinking bracket.\n')
                        end
                        obj.chain(end) = [];
                        scount(n) = scount(n) + 1;
                        if teta < 0
                            tetamin = teta;
                        else
                            tetamax = teta;
                        end
                    end
                    teta = tetamin + rand*(tetamax - tetamin);
                end
                if ~success
                    warning('Possibly degenerate slice.')
                end
                if verb > 0
                    fprintf('Drawing sample %d (of %d)...done. ', n, nsamples)
                    fprintf('Sample scored L = %g.\n', Lfp)
                end
                if ~mod(n,nbackup)
                    obj.dump_chain();
                    %TODO: maybe call reduce_tofs?
                end
            end
            
            warning on TOFPIPE:realize_model_fail
            ET = toc(t_smpl);
            if verb > 0
                fprintf('\nAll samples drawn (%f hit rate). Run time was %s\n',...
                    1/mean(scount), lower(seconds2human(ET)))
            end
            
            % Return average slice hit rate
            hitrate = 1/mean(scount);
        end
    end % End of graveyard block
    
    %% Access methods
    methods
        function set.opts(obj,val)
            obj.opts = pipeset(val);
        end
        
        function set.observables(obj,val)
            if isempty(val), obj.observables = []; return, end
            failmsg = ['Expected scalar struct with fields: ',...
                'a0, M, m, J2, J4'];
            assert(isstruct(val) && isscalar(val), failmsg)
            C = {'a0','M','m','J2','J4'};
            assert(all(isfield(val, C)), failmsg)
            obj.observables = val;
        end
        
        function set.lossfunction(obj,val)
            if isempty(val), obj.lossfunction = []; return, end
            validateattributes(val,{'function_handle'},{'scalar'})
            obj.lossfunction = val;
        end
        
        function set.likefunction(obj,val)
            if isempty(val), obj.likefunction = []; return, end
            validateattributes(val,{'function_handle'},{'scalar'})
            obj.likefunction = val;
        end
        
        function val = get.nlinks(obj)
            val = length(obj.chain);
        end
        
        function val = get.multlinks(obj)
            val = sum([obj.chain.multiplicity]);
        end
        
        function val = get.badlinks(obj)
            if isempty(obj.chain)
                val = 0;
            else
                ind = cellfun(@(x)isempty(x)||isnan(x), {obj.chain.runtime});
                val = sum(ind);
            end
        end
    end % End of access methods block
    
    %% Static methods
    methods (Static)
    end % End of static methods block
end

%% Class-related functions and misc one-liners
function tof = realize_model(mdl, x, mdlopts)
try
    tof = mdl(x);
    if ~isempty(mdlopts)
        tof.opts = mdlopts;
    end
    tof.opts.verbosity = tof.opts.verbosity - 1; % make sense really...
catch ME
    warning('TOFPIPE:realize_model_fail',...
        'Model construction failed with message from %s:\n"%s"',...
        ME.stack(1).name, ME.message)
    tof = [];
end
end

function qr = quick_reject_model(tof)
% Quick reject is yes/no, the value of qr is not meaningful.
if any(tof.rhoi < 0)
    qr = 1;
else
    deltas = diff(tof.rhoi);
    qr = any(deltas < 0);
end
end

function ET = converge_model(tof, popts)
try
    ET = tof.relax_to_HE;
    if popts.matchreq
        tof.si = tof.si*tof.radius/tof.a0;
    end
    if popts.matchmass
        tof.rhoi = tof.rhoi*tof.mass/tof.M;
    end
catch ME
    warning('Model convergence failed with message from %s:\n"%s"',...
        ME.stack(1).name, ME.message)
    ET = NaN;
end
end

function p = basic_parser()
% Return inputParser with common parameters.

p = inputParser();
p.addParameter('axes',[],@(x)isempty(x)||(isscalar(x)&&isgraphics(x, 'axes')))
p.addParameter('cmap',@parula,@(x)isa(x,'function_handle'))
p.addParameter('target','screen',@(x)ischar(x)&&isrow(x))
p.addParameter('lstyle','-')
p.addParameter('color',[])
end

function [fh, ah] = get_canvas(target)
%GET_CANVAS Return handles to new figure and axes set up the way I like.
% This is a copy of the ngraf package function. Try help ngraf.get_canvas for full
% documentation.

% Inputs
narginchk(0,1)
if nargin == 0, target = 'screen'; end
target = validatestring(target, {'screen','projector','publication'});

% For all targets
fh = figure;
set(fh, 'defaultTextInterpreter', 'latex')
set(fh, 'defaultLegendInterpreter', 'latex')
set(fh, 'defaultLegendFontSize', 12)
set(fh, 'defaultLegendFontSizeMode', 'manual')
set(fh, 'defaultLineLinewidth', 2)
ah = axes;
ah.Box = 'on';
hold(ah, 'on');
ah.XLabel.Interpreter = 'latex';
ah.XLabel.FontSize = 14;
ah.YLabel.Interpreter = 'latex';
ah.YLabel.FontSize = 14;
ah.Title.FontSize = 14;

% Center and scale (this often auto-tweaks other stuff)
if ismember(target, {'projector','publication'})
    fh.Position(3:4) = [768, 576];
    ngraf.center_fig(fh);
end

% The grey background is a signature MATLAB visual, but jarring for outsiders
fh.Color = 'w';

% Magnify tick labels (they inherit from axes font)
if ismember(target, {'projector','publication'})
    ah.FontSize = 20;
end

% Magnify axes labels (must come AFTER parent axes font size)
if ismember(target, {'projector','publication'})
    ah.XAxis.Label.FontSize = 24;
    ah.YAxis.Label.FontSize = 26;
end

% Magnify ticks
if ismember(target, {'projector','publication'})
    ah.XAxis.TickLength(1) = 0.02;
    ah.YAxis.TickLength(1) = 0.02;
end

% Thicker axes lines are easier to see, especially in two-column typesetting
if ismember(target, {'publication'})
    ah.XAxis.LineWidth = 2;
    ah.YAxis.LineWidth = 2;
end

% And BYU
end
