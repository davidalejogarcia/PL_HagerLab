function G = Gauss_1Cmp_fun(k0,X,n_total)

mu = k0(1);
sigma = k0(2);
g = normpdf(X,mu,sigma);
g = g./sum(g);
G = g*n_total;


