function [G_Coef, G_Sigma, G_Fit] = Gauss_1Cmp_fit(data, k0)

% Calculate best exponential fit with a the lsqnonlin non linear least
% square fitting routine


% Read input data

tlist =data(:,1);
Y = data(:,2);
% A0 = max(Y);


% define initial values for the parameters
COEF0 = k0;
n_total = sum(Y);

% define lower boundaries for the fitted parameters
COEF_LB = [-Inf, 0];
% define upper boundaries for the fitted parameters
COEF_UB = [Inf, Inf];

% Define anonymous function for the fitting
fitfun = @(COEF) (Gauss_1Cmp_fun(COEF, tlist,n_total) -Y);

% Select Options for the fitting
options = optimset('FunValCheck','off','Display','off');

% run fitting routibe
[G_Coef, resNorm, residuals,exitflag,output,lambda,jacobian] = lsqnonlin(fitfun,COEF0,COEF_LB,COEF_UB,options);
ci = nlparci(G_Coef,residuals,'jacobian',jacobian);
G_Sigma = (ci(:,2) - ci(:,1))/2;
G_Sigma = G_Sigma';

% Compute output
G_Fit(:,1) = data(:,1);
G_Fit(:,2)= Gauss_1Cmp_fun(G_Coef,G_Fit(:,1),n_total);
