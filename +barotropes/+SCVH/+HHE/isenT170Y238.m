function eos = isenT170Y238()
%ISENT170Y238 SCVH H/He adiabat; T_1bar=170 K; Y=0.238.

tblfile = [mfilename('fullpath'),'.adiabat'];
eos = barotropes.SCVH.SCVH(tblfile);
eos.T_1bar = 170;
eos.Y = 0.238;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
