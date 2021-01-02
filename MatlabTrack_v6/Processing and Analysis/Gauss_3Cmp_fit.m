function [G_Coef, G_Sigma, G_Fit] = Gauss_3Cmp_fit(data, k0)

% Calculate best exponential fit with a the lsqnonlin non linear least
% square fitting routine


% Read input data

tlist =data(:,1);
Y = data(:,2);
f01 = 0.5;
f02 = 0.25;
% A0 = max(Y);


% define initial values for the parameters
COEF0 = [k0,f01,f02];
n_total = sum(Y);

% define lower boundaries for the fitted parameters
COEF_LB = [-Inf, 0, -Inf, 0, -Inf, 0, 0, 0];
% define upper boundaries for the fitted parameters
COEF_UB = [Inf, Inf, Inf, Inf, Inf, Inf, 1, 1];

% Define anonymous function for the fitting
fitfun = @(COEF) (Gauss_3Cmp_fun(COEF, tlist,n_total) -Y);

% Select Options for the fitting
options = optimset('FunValCheck','off','Display','off');

% run fitting routibe
[G_Coef, resNorm, residuals,exitflag,output,lambda,jacobian] = lsqnonlin(fitfun,COEF0,COEF_LB,COEF_UB,options);
ci = nlparci(G_Coef,residuals,'jacobian',jacobian);
G_Sigma = (ci(:,2) - ci(:,1))/2;
G_Sigma = G_Sigma';

% Compute output
G_Fit(:,1) = data(:,1);
G_Fit(:,2)= Gauss_3Cmp_fun(G_Coef,G_Fit(:,1),n_total);
G_Fit(:,3) = Gauss_1Cmp_fun(G_Coef(1:2),G_Fit(:,1),G_Coef(7)*n_total);
G_Fit(:,4) = Gauss_1Cmp_fun(G_Coef(3:4),G_Fit(:,1),G_Coef(8)*n_total);
G_Fit(:,5) = Gauss_1Cmp_fun(G_Coef(5:6),G_Fit(:,1),(1 - G_Coef(7) - G_Coef(8))*n_total);


