function [Esp_Coef, Esp_Sigma, Esp_Fit] = BleachingAndRecovery_fit(data,k0s,val0)

% Calculate best exponential fit with a the lsqnonlin non linear least
% square fitting routine


% Read input data

tlist =data(:,1);
Y = data(:,2);
A0 = max(Y);


% define initial values for the parameters
COEF0 = [k0s(1),k0s(2), k0s(3), A0];

% define lower boundaries for the fitted parameters
COEF_LB = [0, 0, 0, 0];
% define upper boundaries for the fitted parameters
COEF_UB = [Inf,Inf, Inf, Inf];

% Define anonymous function for the fitting
fitfun = @(COEF) (BleachingAndRecovery_fun(COEF, tlist,val0) -Y);

% Select Options for the fitting
options = optimset('FunValCheck','off','Display','off');

% run fitting routibe
[Esp_Coef, resNorm, residuals,exitflag,output,lambda,jacobian] = lsqnonlin(fitfun,COEF0,COEF_LB,COEF_UB,options);
ci = nlparci(Esp_Coef,residuals,'jacobian',jacobian);
Esp_Sigma = (ci(:,2) - ci(:,1))/2;
Esp_Sigma = Esp_Sigma';

% Compute output
Esp_Fit(:,1) = data(:,1);
Esp_Fit(:,2)= BleachingAndRecovery_fun(Esp_Coef,Esp_Fit(:,1),val0);

