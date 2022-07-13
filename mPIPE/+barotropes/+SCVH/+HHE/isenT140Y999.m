function eos = isenT140Y999()
%ISENT140Y999 SCVH H/He adiabat; T_1bar=140 K; Y=1.0 (pure He).

tblfile = [mfilename('fullpath'),'.adiabat'];
eos = barotropes.SCVH.SCVH(tblfile);
eos.T_1bar = 140;
eos.Y = 1.0;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
