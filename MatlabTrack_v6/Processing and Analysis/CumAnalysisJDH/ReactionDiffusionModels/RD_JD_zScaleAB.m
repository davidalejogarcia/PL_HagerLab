function zScale = RD_JD_zScaleAB(deltaZ, PAR, flag)

if flag == 1;
tlist = PAR(1,1);
else
tlist = PAR(:,1);
end
D = PAR(1,2);
n = PAR(2,2);

COEF(1) = D;
COEF(2) = deltaZ;



IntFun = @(z)RD_JD_zScaleAB_int(z, tlist, COEF);
LowLimit = - deltaZ/2;
HighLimit = deltaZ/2;
zScale=n*quadv(IntFun, LowLimit, HighLimit)/deltaZ;

