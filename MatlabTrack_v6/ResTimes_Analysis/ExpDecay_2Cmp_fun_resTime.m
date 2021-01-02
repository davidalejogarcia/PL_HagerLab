function esp = ExpDecay_2Cmp_fun_resTime(coeff, x)


k1 = coeff(1);
k2 = coeff(2);
f = coeff(3);
Ceq = coeff(4);


esp = Ceq*(f*k1*(exp(-k1*x)) + (1-f)*k2*exp(-k2*x));