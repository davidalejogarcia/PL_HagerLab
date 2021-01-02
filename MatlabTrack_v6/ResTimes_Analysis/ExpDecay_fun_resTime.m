function esp = ExpDecay_fun_resTime(coeff, x)


k = coeff(1);
a= coeff(2);



esp = a*k*(exp(-k*x));


