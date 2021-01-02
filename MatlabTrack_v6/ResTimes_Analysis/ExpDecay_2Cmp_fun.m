function esp = ExpDecay_2Cmp_fun(coeff, x)


k1 = coeff(1);
k2 = coeff(2);
f = coeff(3);
a = coeff(4);


esp = a*(f*(exp(-k1*x)) + (1-f)*exp(-k2*x));