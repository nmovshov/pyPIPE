function eos = isenT140Y180()
%ISENT140Y180 REOS3 H/He adiabat; T_1bar=140 K; Y=0.180.

tblfile = [mfilename('fullpath'),'.dat'];
eos = barotropes.REOS3.REOS3(tblfile);
eos.T_1bar = 140;
eos.Y = 0.180;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
