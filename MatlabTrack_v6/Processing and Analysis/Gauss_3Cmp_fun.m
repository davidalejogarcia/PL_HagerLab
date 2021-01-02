function G = Gauss_3Cmp_fun(k0,X,n_total)

mu1 = k0(1);
sigma1 = k0(2);
mu2 = k0(3);
sigma2 = k0(4);
mu3 = k0(5);
sigma3 = k0(6);

frac1 = k0(7);
frac2 = k0(8);


g1 = normpdf(X,mu1,sigma1);
g2 = normpdf(X,mu2,sigma2);
g3 = normpdf(X,mu3,sigma3);

g = frac1*g1 + frac2*g2 + (1-frac1-frac2)*g3;
g = g./sum(g);
G = n_total*g;
