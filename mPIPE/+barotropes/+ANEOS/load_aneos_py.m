function eos = load_aneos_py(material)
%LOAD_ANEOS_PY Return handle to Chris's python aneos.eos class.

if nargin < 1 || isempty(material), material = 'serpentine'; end
[pp,~,~] = fileparts(mfilename('fullpath'));
eos = py.aneos.eos(pp, material);
end
