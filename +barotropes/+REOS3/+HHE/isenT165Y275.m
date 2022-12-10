function eos = isenT165Y275()
%ISENT165Y275 REOS3 H/He adiabat; T_1bar=165 K; Y=0.275.

tblfile = [mfilename('fullpath'),'.dat'];
eos = barotropes.REOS3.REOS3(tblfile);
eos.T_1bar = 165;
eos.Y = 0.275;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
