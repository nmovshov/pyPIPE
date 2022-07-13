function eos = isenP10T150Y275()
%ISENP10T140Y275 SCVH H/He adiabat; T_10bar=150 K; Y=0.275.

tblfile = [mfilename('fullpath'),'.adiabat'];
eos = barotropes.SCVH.SCVH(tblfile);
eos.T_1bar = 70;
eos.Y = 0.275;
eos.X = 1-eos.Y;
eos.name = mfilename;
end
