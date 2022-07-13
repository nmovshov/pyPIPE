function eos = serpentine_on_isenT140Y275()
%SERPENTINE_ON_ISENT140Y275 ANEOS serpentine density on SCVH H/He adiabat.

tblfile = [mfilename('fullpath'),'.adiabat'];
eos = barotropes.SCVH.SCVH(tblfile);
eos.T_1bar = 140;
eos.Y = 0.275;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
