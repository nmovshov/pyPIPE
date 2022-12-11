function eos = isenT140Y125()
%ISENT140Y012 SCVH H/He adiabat; T_1bar=140 K; Y=0.125.

tblfile = [mfilename('fullpath'),'.adiabat'];
eos = barotropes.SCVH.SCVH(tblfile);
eos.T_1bar = 140;
eos.Y = 0.125;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
