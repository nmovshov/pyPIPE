function eos = serpentine_on_isenP10T150Y275()
%SERPENTINE_ON_ISENP10T150Y275 ANEOS serpentine density on SCVH H/He adiabat.

tblfile = [mfilename('fullpath'),'.adiabat'];
eos = barotropes.SCVH.SCVH(tblfile);
eos.T_1bar = 70;
eos.Y = 0.275;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
