function zh = zoom_inset(ah, zoomed_xlim, zoomed_ylim)
%ZOOM_INSET Place a copy of an axes in the top right corner of parent figure.
%   ZOOM_INSET(ah) makes a copy of the axes ah, removes the axis labels,
%   legend, and title, reduces its size and places it flush with the top right
%   corner of the parent figure. Assuming that the parent figure had one axes
%   to begin with the result is a well-placed "inset" useful, e.g., for
%   emphasizing a zoomed-in region of the plotted range. A handle to the new
%   copy is returned; there will almost always be a need to further customize
%   aspects of the inset axes.
%
%   ZOOM_INSET(ah, zoomed_xlim, zoomed_ylim) also updates the inset's axis limits.
%
% Example:
%   plot(-pi:0.1:pi, sin(-pi:0.1:pi))
%   zoom_inset(gca, [pi/2-1, pi/2+1])

narginchk(0,3)
if nargin == 0
    fprintf('Usage:\n\tzoom_inset(ah,zoomed_xlim,zoomed_ylim).\n')
    return
end
if nargin < 2 , zoomed_xlim = []; end
if nargin < 3 , zoomed_ylim = []; end
validateattributes(ah,{'matlab.graphics.axis.Axes'},{'scalar'})

zh = copyobj(ah, ah.Parent);
zh.XLabel.reset;
zh.YLabel.reset;
zh.Title.reset;
zh.Legend.reset;
zh.Position(3:4) = 0.6*ah.Position(3:4);
zh.Position(1:2) = ah.Position(1:2) + 0.4*ah.Position(3:4);
if ~isempty(zoomed_xlim), zh.XLim = zoomed_xlim; end
if ~isempty(zoomed_ylim), zh.YLim = zoomed_ylim; end

end
