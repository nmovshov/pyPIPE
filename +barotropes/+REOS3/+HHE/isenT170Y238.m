function eos = isenT170Y238()
%ISENT170Y238 REOS3 H/He adiabat; T_1bar=170 K; Y=0.238.

tblfile = [mfilename('fullpath'),'.dat'];
eos = barotropes.REOS3.REOS3(tblfile);
eos.T_1bar = 170;
eos.Y = 0.238;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
