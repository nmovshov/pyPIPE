function [fh, ah] = polyprofs_grid(filename,varargin)
% Comparison grid for polynomial models with varying degree and precision.

%% Inputs
if nargin == 0
    fprintf('Usage:\n\tpolyprofs_grid(filename,''name'',''value'')\n')
    return
end
narginchk(1,inf)
p = plots.basic_parser();
p.parse(varargin{:})
pr = p.Results;

%% Get a big canvas for this multi-axes figure
[fh, ah] = ngraf.get_canvas(pr.target);
delete(ah)
fh.WindowState = 'max';

%% Load the data
data = load(filename);

%% I'll go through the grid manually, for now
ah = gobjects(1,8);

subplot(4,2,1)
plots.average_of_profs(data.d5x1, 'ax', gca, 'leg', false);
textbp('d5x1')

subplot(4,2,2)
plots.average_of_profs(data.d5x10, 'ax', gca, 'leg', false);
textbp('d5x10')

subplot(4,2,3)
plots.average_of_profs(data.d6x1, 'ax', gca, 'leg', false);
textbp('d6x1')

subplot(4,2,4)
plots.average_of_profs(data.d6x10, 'ax', gca, 'leg', false);
textbp('d6x10')

subplot(4,2,5)
plots.average_of_profs(data.d7x1, 'ax', gca, 'leg', false);
textbp('d7x1')

subplot(4,2,6)
plots.average_of_profs(data.d7x10, 'ax', gca, 'leg', false);
textbp('d7x10')

subplot(4,2,7)
plots.average_of_profs(data.d8x1, 'ax', gca, 'leg', false);
textbp('d8x1')

subplot(4,2,8)
plots.average_of_profs(data.d8x10, 'ax', gca, 'leg', false);
textbp('d8x10')

%% BYU
end
