function eos = He_on_isenT140Y111()
%HE_ON_ISENT140Y111 SCVH Helium density on SCVH H/He adiabat.

tblfile = [mfilename('fullpath'),'.adiabat'];
eos = barotropes.SCVH.SCVH(tblfile);
eos.T_1bar = 140;
eos.Y = 0.111;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
