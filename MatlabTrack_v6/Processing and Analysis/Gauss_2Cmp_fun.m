function G = Gauss_2Cmp_fun(k0,X,n_total)

mu1 = k0(1);
sigma1 = k0(2);
mu2 = k0(3);
sigma2 = k0(4);
frac1 = k0(5);


g1 = normpdf(X,mu1,sigma1);
g2 = normpdf(X,mu2,sigma2);

g = frac1*g1 + (1-frac1)*g2;
g = g./sum(g);
G = n_total*g;
