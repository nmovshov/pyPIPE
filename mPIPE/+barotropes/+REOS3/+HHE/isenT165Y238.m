function eos = isenT165Y238()
%ISENT165Y238 REOS3 H/He adiabat; T_1bar=165 K; Y=0.238.

tblfile = [mfilename('fullpath'),'.dat'];
eos = barotropes.REOS3.REOS3(tblfile);
eos.T_1bar = 165;
eos.Y = 0.238;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
