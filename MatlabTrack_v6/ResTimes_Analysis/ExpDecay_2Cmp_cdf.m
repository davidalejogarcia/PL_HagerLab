function esp = ExpDecay_2Cmp_cdf(x,k1,k2,f)


% k1 = coeff(1);
% k2 = coeff(2);
% f = coeff(3);
% a = coeff(4);


% esp = 1 - (f*(exp(-k1*x)) + (1-f)*exp(-k2*x));
esp = f*expcdf(x,1/k1) + (1-f)*expcdf(x,1/k2);