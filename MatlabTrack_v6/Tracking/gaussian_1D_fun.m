
function I = gaussian_1D(COEF, rlist)

I0=COEF(1);
sigma = COEF(2);
r0= COEF(3);
bkg = COEF(4);

r = rlist;

I=I0.*exp(-(r-r0).^2/(2*sigma^2))+bkg;