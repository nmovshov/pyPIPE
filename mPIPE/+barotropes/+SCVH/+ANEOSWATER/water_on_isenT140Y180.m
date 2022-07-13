function eos = water_on_isenT140Y180()
%WATER_ON_ISENT140Y180 ANEOS water density on SCVH H/He adiabat.

tblfile = [mfilename('fullpath'),'.adiabat'];
eos = barotropes.SCVH.SCVH(tblfile);
eos.T_1bar = 140;
eos.Y = 0.180;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
