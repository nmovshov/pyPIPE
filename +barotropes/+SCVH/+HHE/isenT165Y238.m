function eos = isenT165Y238()
%ISENT165Y238 SCVH H/He adiabat; T_1bar=165 K; Y=0.238.

tblfile = [mfilename('fullpath'),'.adiabat'];
eos = barotropes.SCVH.SCVH(tblfile);
eos.T_1bar = 165;
eos.Y = 0.238;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
