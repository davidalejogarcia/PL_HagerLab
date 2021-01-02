function Integrand = RD_JD_zScaleAB_int(z, t, COEF)


% Read input

D = COEF(1);
deltaZ = COEF(2);

n = 40;
n = 0: n;


[n, t] =  meshgrid(n,t);


psi = sqrt(4*D*t);

Integrand = (-1).^n.*(erfc((z + (2*n + 1)*deltaZ/2)./psi) + ...
    erfc((- z + (2*n + 1)*deltaZ/2)./psi));

Integrand = 1 - sum(Integrand, 2);