function [fh, ah] = get_canvas(target)
%GET_CANVAS Return handles to new figure and axes set up the way I like.
%   [fh, ah] = GET_CANVAS() creates and returns handles to a new figure and
%   axes and modifies the default values of a few properties to get closer to
%   publication/presentation style. Typically the result will have to be
%   tweaked a little more for the actual camera-ready figure, but in ways that
%   will depend on the content and can't be generalized. On the other hand, I
%   don't want to override the groot defaults even of properties I virtually
%   always modify. This is partly because I find that the defaults work well
%   for the quick glance, exploration phase, and partly because I want to stay
%   aware of what the defaults are.
%
%   [fh, ah] = GET_CANVAS(target) further modifies the defaults based on where
%   the figure is destined to end up, as defined by the string target.
%   Currently the supported targets are 'screen', 'publication', and
%   'projector', and in fact 'publication' and 'projector' have almost the same
%   look, but this can change. The default is target='screen'. (Hint: partial
%   match works.)

%% Inputs
narginchk(0,1)
if nargin == 0, target = 'screen'; end
target = validatestring(target, {'screen','projector','publication'});

%% For all targets
fh = figure;
set(fh, 'defaultTextInterpreter', 'latex')
set(fh, 'defaultLegendInterpreter', 'latex')
set(fh, 'defaultLineLinewidth', 2)
ah = axes;
ah.Box = 'on';
hold(ah, 'on');
ah.XLabel.Interpreter = 'latex';
ah.XLabel.FontSize = 14;
ah.YLabel.Interpreter = 'latex';
ah.YLabel.FontSize = 14;
ah.Title.FontSize = 14;

%% Center and scale (this often auto-tweaks other stuff)
if ismember(target, {'projector','publication'})
    fh.Position(3:4) = [768, 576];
    ngraf.center_fig(fh);
end

%% The grey background is a signature MATLAB visual, but jarring for outsiders
fh.Color = 'w';

%% Magnify tick labels (they inherit from axes font)
if ismember(target, {'projector','publication'})
    ah.FontSize = 20;
end

%% Magnify axes labels (must come AFTER parent axes font size)
if ismember(target, {'projector','publication'})
    ah.XAxis.Label.FontSize = 24;
    ah.YAxis.Label.FontSize = 26;
end

%% Magnify ticks
if ismember(target, {'projector','publication'})
    ah.XAxis.TickLength(1) = 0.02;
    ah.YAxis.TickLength(1) = 0.02;
end

%% Thicker axes lines are easier to see, especially in two-column typesetting
if ismember(target, {'publication'})
    ah.XAxis.LineWidth = 2;
    ah.YAxis.LineWidth = 2;
end

%% And BYU
return

end
