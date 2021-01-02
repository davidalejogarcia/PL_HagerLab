function esp = GeoDecay_fun(coeff, x)


k = coeff(1);
a= coeff(2);



% esp = a*(exp(-k*x));

esp = a*(1-k).^x;
