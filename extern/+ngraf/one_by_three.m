function fh = one_by_three(fh)
%ONE_BY_THREE Adjust position of child axes to 1-by-3 grid.
%   ONE_BY_THREE() repositions the children of the current figure to a 1-by-3
%   grid. Of course, there should be exactly 3 child axes.
%
%   ONE_BY_THREE(fh) positions the children of fh.

narginchk(0,1)
if nargin == 0, fh = gcf; end
validateattributes(fh,{'matlab.ui.Figure'},{'scalar'})

ah = flipud(findobj(fh, 'type', 'axes'));
assert(length(ah) == 3, '1-by-3 grid has 3 spaces.')

ah(1).Position = [0.1300, 0.1100, 0.2134, 0.8150];
ah(2).Position = [0.4108, 0.1100, 0.2134, 0.8150];
ah(3).Position = [0.6916, 0.1100, 0.2134, 0.8150];
end
