function fh = one_by_two(fh)
%ONE_BY_TWO Adjust position of child axes to 1-by-2 grid.
%   ONE_BY_TWO() repositions the children of the current figure to a 1-by-2
%   grid. Of course, there should be exactly 2 child axes.
%
%   ONE_BY_TWO(fh) positions the children of fh.

narginchk(0,1)
if nargin == 0, fh = gcf; end
validateattributes(fh,{'matlab.ui.Figure'},{'scalar'})

ah = flipud(findobj(fh, 'type', 'axes'));
assert(length(ah) == 2, '1-by-2 grid has 2 spaces.')

ah(1).Position = [0.1300, 0.1100, 0.3347, 0.8150];
ah(2).Position = [0.5703, 0.1100, 0.3347, 0.8150];
end
