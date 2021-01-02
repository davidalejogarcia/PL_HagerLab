function esp = ExpDecay_3Cmp_fun(coeff, x)


k1 = coeff(1);
k2 = coeff(2);
k3 = coeff(3);
f1 = coeff(4);
f2 = coeff(5);
a = coeff(6);



esp = a*(f1*(exp(-k1*x)) + f2*exp(-k2*x) + (1 - f1 - f2)*exp(-k3*x));