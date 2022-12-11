function eos = isenT140Y000()
%ISENT140Y000 SCVH H/He adiabat; T_1bar=140 K; Y=0.0 (pure H).

tblfile = [mfilename('fullpath'),'.adiabat'];
eos = barotropes.SCVH.SCVH(tblfile);
eos.T_1bar = 140;
eos.Y = 0.0;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
