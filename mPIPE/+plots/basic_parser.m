function p = basic_parser()
% Return inputParser with common parameters.

p = inputParser();
p.addParameter('axes',[],@(x)isempty(x)||(isscalar(x)&&isgraphics(x, 'axes')))
p.addParameter('cmap',@parula,@(x)isa(x,'function_handle'))
p.addParameter('target','screen',@(x)ischar(x)&&isrow(x))
p.addParameter('lstyle','-')
p.addParameter('color',[])
p.addParameter('legend',true,@(x)isscalar(x)&&islogical(x))
p.addParameter('title','',@(x)ischar(x)&&isrow(x))
p.addParameter('label','',@(x)ischar(x)&&isrow(x))
end
