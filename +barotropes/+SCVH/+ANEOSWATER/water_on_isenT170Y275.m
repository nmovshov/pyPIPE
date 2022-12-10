function eos = water_on_isenT170Y275()
%WATER_ON_ISENT170Y275 ANEOS water density on SCVH H/He adiabat.

tblfile = [mfilename('fullpath'),'.adiabat'];
eos = barotropes.SCVH.SCVH(tblfile);
eos.T_1bar = 170;
eos.Y = 0.275;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
