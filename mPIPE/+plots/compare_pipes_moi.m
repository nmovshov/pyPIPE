function compare_pipes_moi(pipes, varargin)
% COMPARE_PIPES_MOI Compare MoI histograms from different samples.

if nargin == 0
    fprintf('Usage:\n\tcompare_pipes_moi(pipes, ''name'',value)\n')
    fprintf('\t\tpipes is an array of TOFPIPEs.\n')
    return
end

validateattributes(pipes,{'TOFPIPE'},{'vector'})
p = inputParser;
p.addParameter('target','proj')
p.addParameter('alfa',0.8)
p.parse(varargin{:})
pr = p.Results;

[~,~] = ngraf.get_canvas(pr.target);
for k=1:length(pipes)
    moi = double(pipes(k).getmdlfield('NMoI'));
    hh = histogram(moi,'Normalization','pdf');
    hh.DisplayName = pipes(k).name;
end
ylabel('P density')
legend
end
