function esp = ExpDecay_fun(coeff, x)


k = coeff(1);
a= coeff(2);



esp = a*(exp(-k*x));


