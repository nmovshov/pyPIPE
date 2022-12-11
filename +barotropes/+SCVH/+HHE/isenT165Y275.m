function eos = isenT165Y275()
%ISENT165Y275 SCVH H/He adiabat; T_1bar=165 K; Y=0.275.

tblfile = [mfilename('fullpath'),'.adiabat'];
eos = barotropes.SCVH.SCVH(tblfile);
eos.T_1bar = 165;
eos.Y = 0.275;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
