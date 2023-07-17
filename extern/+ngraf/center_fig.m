function fh = center_fig(fh)
%CENTER_FIG Align center of figure with center of screen.
%   CENTER_FIG() positions the current figure at the center of the main display.
%
%   CENTER_FIG(fh) positions the figure identified by handle fh.

narginchk(0,1)
if nargin == 0, fh = gcf; end
validateattributes(fh,{'matlab.ui.Figure'},{'scalar'})

gh = groot;
fh.Position(1) = round(gh.ScreenSize(3)/2) - round(fh.Position(3)/2);
fh.Position(2) = round(gh.ScreenSize(4)/2) - round(fh.Position(4)/2);
end
