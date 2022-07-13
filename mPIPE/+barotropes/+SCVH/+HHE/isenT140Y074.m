function eos = isenT140Y074()
%ISENT140Y074 SCVH H/He adiabat; T_1bar=140 K; Y=0.074.

tblfile = [mfilename('fullpath'),'.adiabat'];
eos = barotropes.SCVH.SCVH(tblfile);
eos.T_1bar = 140;
eos.Y = 0.074;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
