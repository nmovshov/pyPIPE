function compare_pipes_profiles(pipes, varargin)
% COMPARE_PIPES_PROFILES Compare rho(s) profiles from different samples.

global proot % project root defined in setws

if nargin == 0
    fprintf('Usage:\n\tcompare_pipes_profiles(pipes, ''name'',value)\n')
    fprintf('\t\tpipes is an array of TOFPIPEs.\n')
    return
end

validateattributes(pipes,{'TOFPIPE'},{'vector'})
p = inputParser;
p.addParameter('savehc',false)
p.addParameter('target','proj')
p.addParameter('titles',true)
p.addParameter('inset',true)
p.addParameter('sbys',true)
p.addParameter('linky',true)
p.addParameter('nlines',30)
p.addParameter('alfa',0.8)
p.parse(varargin{:})
pr = p.Results;

if pr.sbys
    assert(length(pipes) <= 3, 'side-by-side compares up to 3 chains.')
    pr.inset = false;
    sbyf = ngraf.get_canvas(pr.target);
    delete(sbyf.Children);
    sbyf.WindowState = 'max';
end

warning('off','MATLAB:handle_graphics:exceptions:SceneNode')
for k=1:length(pipes)
    pipe = pipes(k);
    prof = pipe.getmdlfield('rhoi');
    svec = pipe.chain(1).tof.si;
    
    [fh,ah,~] = plots.ensemble_of_profs(prof,'svec',svec,'tar',pr.target,...
        'nlines',pr.nlines,'alfa',pr.alfa);
    legend off
    if pr.titles
        title(pipe.name)
    end
    
    if pr.sbys, copyobj(ah,sbyf); end
    
    if pr.inset
        [tfh,tah,~] = plots.ensemble_of_profs(prof,'svec',svec,'tar',pr.target,...
            'nlines',pr.nlines,'alfa',pr.alfa);
        tah.XLim(1) = 0.6;
        tah.YLimMode = 'auto';
        tah.XTick = 0.6:0.1:0.9;
        tah.YMinorTick = 'on';
        tah.XLabel.reset;
        tah.YLabel.reset;
        pause(1)
        tah = copyobj(tah, fh);
        close(tfh)
        tah.Position(3:4) = 0.6*ah.Position(3:4);
        tah.Position(1:2) = ah.Position(1:2) + 0.4*ah.Position(3:4);
    end
    
    if pr.sbys, close(fh); end
end

if pr.sbys
    if length(pipes) == 2, ngraf.one_by_two(sbyf); end
    if length(pipes) == 3, ngraf.one_by_three(sbyf); end
    if pr.linky
        linkaxes(sbyf.Children);
    end
end
end
