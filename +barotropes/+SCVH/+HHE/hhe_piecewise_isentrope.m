function eos = hhe_piecewise_isentrope(yvec, pvec, t1vec)
%HHE_PIECEWISE_ISENTROPE Multiple SCVH H/He stitched adiabats.

if nargin == 0
    fprintf('Usage:\n\thhe_piecewise_isentrope(yvec, pvec, t1vec)\n')
    return
end
narginchk(3,3)
validateattributes(yvec,{'numeric'},{'row','>=',0,'<=',1})
validateattributes(pvec,{'numeric'},{'row','>=',0})
validateattributes(t1vec,{'numeric'},{'row','>=',0})
assert(isequal(size(yvec),size(pvec),size(t1vec)))

%TODO: implement

eos = barotropes.SCVH.SCVH();
eos.T_1bar = t1vec(1);
eos.Y = yvec;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
