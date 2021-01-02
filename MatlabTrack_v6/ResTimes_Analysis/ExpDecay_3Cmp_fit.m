function [Esp_Coef, Esp_Sigma Esp_Fit] = ExpDecay_3Cmp_fit(data, k0s)

% Calculate best two compoonent exponential fit with a the lsqnonlin non 
%linear least square fitting routine.


% Read input data

tlist =data(:,1);
Y = data(:,2);


k01 = k0s(1);
k02 = k0s(2);
k03 = k0s(3);
f01 = 0.5;
f02 = 0.25;
A0 = max(Y);


% define initial values for the parameters
COEF0 = [k01, k02, k03, f01, f02, A0];


% define lower boundaries for the fitted parameters
COEF_LB = [0, 0, 0 , 0, 0, 0];
%COEF_LB = [0 , 0];
% define upper boundaries for the fitted parameters
COEF_UB = [Inf, Inf, Inf,1, 1,Inf];
% COEF_UB = [1,Inf];
% Define anonymous function for the fitting
fitfun = @(COEF) (ExpDecay_3Cmp_fun(COEF, tlist)) - (Y);

% Select Options for the fitting
options = optimset('FunValCheck','off', 'MaxIter', 1000, 'MaxFunEvals', 1000);

% run fitting routibe
[Esp_Coef, resNorm, residuals,exitflag,output,lambda,jacobian] = lsqnonlin(fitfun,COEF0,COEF_LB,COEF_UB,options);


ci = nlparci(Esp_Coef,residuals,'jacobian',jacobian);
Esp_Sigma = (ci(:,2) - ci(:,1))/2;
Esp_Sigma = Esp_Sigma';



% Compute output
Esp_Fit(:,1) = data(:,1);
Esp_Fit(:,2) = ExpDecay_3Cmp_fun(Esp_Coef, data(:,1));

%Compute the individual components --DB
% Esp_Coef_1 = [Esp_Coef(1) Esp_Coef(4)*Esp_Coef(3)];
Esp_Coef_1 = [Esp_Coef(1) Esp_Coef(4)*Esp_Coef(6)];
Esp_Fit(:,3) = ExpDecay_fun(Esp_Coef_1,data(:,1));

% Esp_Coef_2 = [Esp_Coef(2) Esp_Coef(4)*(1 - Esp_Coef(3))];
Esp_Coef_2 = [Esp_Coef(2) Esp_Coef(5)*Esp_Coef(6)];
Esp_Fit(:,4) = ExpDecay_fun(Esp_Coef_2,data(:,1));

Esp_Coef_3 = [Esp_Coef(3) Esp_Coef(6)*(1 - Esp_Coef(4) - Esp_Coef(5))];
Esp_Fit(:,5) = ExpDecay_fun(Esp_Coef_3,data(:,1));