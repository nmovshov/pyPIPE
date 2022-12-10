function eos = isenT140Y240()
%ISENT140Y240 SCVH H/He adiabat; T_1bar=140 K; Y=0.240.

tblfile = [mfilename('fullpath'),'.adiabat'];
eos = barotropes.SCVH.SCVH(tblfile);
eos.T_1bar = 140;
eos.Y = 0.240;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
