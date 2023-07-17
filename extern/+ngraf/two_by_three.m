function fh = two_by_three(fh,mxit)
%TWO_BY_THREE Adjust position of child axes to 2-by-3 grid.
%   TWO_BY_THREE() repositions the children of the current figure to a 2-by-3
%   grid. Of course, there should be exactly 6 child axes.
%
%   TWO_BY_THREE(fh) positions the children of fh.
%
%   TWO_BY_THREE(fh,true) first maximizes the figure window.

narginchk(0,2)
if nargin == 0, fh = gcf; end
if nargin < 2, mxit = true; end
validateattributes(fh,{'matlab.ui.Figure'},{'scalar'})
validateattributes(mxit,{'logical'},{'scalar'})

if mxit, fh.WindowState = 'maximized'; end
    
ah = flipud(findobj(fh, 'type', 'axes'));
assert(length(ah) == 6, '2-by-3 grid has 6 spaces.')

ah(1).Position = [0.1300, 0.5838, 0.2134, 0.3412];
ah(2).Position = [0.4108, 0.5838, 0.2134, 0.3412];
ah(3).Position = [0.6916, 0.5838, 0.2134, 0.3412];
ah(4).Position = [0.1310, 0.1100, 0.2134, 0.3412];
ah(5).Position = [0.4108, 0.1100, 0.2134, 0.3412];
ah(6).Position = [0.6916, 0.1100, 0.2134, 0.3412];
end
