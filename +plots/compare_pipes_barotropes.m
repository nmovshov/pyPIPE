function compare_pipes_barotropes(pipes, varargin)
% COMPARE_PIPES_BAROTROPE Compare rho(P) envelopes from different samples.

global proot % project root defined in setws

if nargin == 0
    fprintf('Usage:\n\tcompare_pipes_barotropes(pipes, ''name'',value)\n')
    fprintf('\t\tpipes is an array of TOFPIPEs.\n')
    return
end

validateattributes(pipes,{'TOFPIPE'},{'vector'})
p = plots.basic_parser;
p.addParameter('savehc',false)
p.addParameter('alfa',1.0,@(x)isscalar(x)&&x>0)
p.addParameter('titles',true)
p.addParameter('inset',true)
p.addParameter('sbys',true)
p.addParameter('linky',true)
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
    pees = pipe.getmdlfield('Pi')/1e11;
    rhos = pipe.getmdlfield('rhoi')/1000;
    
    [fh,ah,~] = plots.envelope_of_barotropes(pees,rhos,target=pr.target,...
        alfa=pr.alfa,color=pr.color);
    legend off
    if pr.titles
        title(pipe.name)
    end
    
    if pr.sbys, copyobj(ah,sbyf); end
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
