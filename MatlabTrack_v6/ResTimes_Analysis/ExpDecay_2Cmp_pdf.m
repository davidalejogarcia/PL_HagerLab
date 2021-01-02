function esp = ExpDecay_2Cmp_pdf(x,k1,k2,f)


% k1 = coeff(1);
% k2 = coeff(2);
% f = coeff(3);
% a = coeff(4);
% x1 = unique(x);
% dx = x1(2) - x1(1);
% 
% esp = (f*(k1*exp(-k1*x)) + (1-f)*exp(-k2*x));
% esp = esp./(sum(esp)*dx);

esp = f*exppdf(x,1/k1) + (1-f)*exppdf(x,1/k2);