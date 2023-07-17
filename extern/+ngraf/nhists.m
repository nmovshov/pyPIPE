function [fh,ah,hh] = nhists(data,data2)
%NHISTS Convenience wrapper to histogram with my favorite defaults.
%   NHISTS(data) where data is a cell array just calls the built-in function
%   histogram successively on data{i}, making sure hold is on and setting
%   normalization to 'pdf', more useful than the default 'count' in my opinion,
%   especially when comparing data sets.
%
%   NHISTS(data1,data2) where data1 and data2 are vectors is an alternative
%   syntax, convenient for comparing two data sets.

% Inputs
if nargin == 0
    fprintf('Usage:\n\tnhists(data1,data2)\n')
    fprintf('where data1 and data2 are vectors, or\n')
    fprintf('\tnhists(data)\n')
    fprintf('where data is a cell array of vectors.\n')
    return
end
narginchk(1,2)
if nargin > 1 % hack for common case of two data sets
    validateattributes(data,{'numeric'},{'vector'},'','data1')
    validateattributes(data2,{'numeric'},{'vector'},'','data2')
    data = {data,data2};
end
if ~iscell(data), data = {data}; end % hack for single dataset input
validateattributes(data,{'cell'},{'vector'},'','data')
assert(all(cellfun(@isvector,data)), 'nhists works on one-dimensional data.')

% Call histogram on each data set
[fh,ah] = ngraf.get_canvas(); % holds ah, fyi
hh = gobjects(size(data));
for k=1:length(data)
    hh(k) = histogram(data{k},'Normalization','pdf');
end

% Some styling and annotation
ylabel('P density')
end
