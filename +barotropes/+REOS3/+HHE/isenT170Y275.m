function eos = isenT170Y275()
%ISENT170Y275 REOS3 H/He adiabat; T_1bar=170 K; Y=0.275.

tblfile = [mfilename('fullpath'),'.dat'];
eos = barotropes.REOS3.REOS3(tblfile);
eos.T_1bar = 170;
eos.Y = 0.275;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
