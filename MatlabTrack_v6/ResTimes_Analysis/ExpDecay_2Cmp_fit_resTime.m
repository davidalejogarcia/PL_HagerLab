function [Esp_Coef, Esp_Sigma Esp_Fit] = ExpDecay_2Cmp_fit_resTime(data, k0s)

% Calculate best two compoonent exponential fit with a the lsqnonlin non 
%linear least square fitting routine.


% Read input data

tlist =data(:,1);
Y = data(:,2);


k01 = k0s(1);
k02 = k0s(2);
f01 = 0.5;
A0 = max(Y)/k01;


% define initial values for the parameters
COEF0 = [k01, k02, f01, A0];


% define lower boundaries for the fitted parameters
COEF_LB = [0, 0, 0 , 0];

% define upper boundaries for the fitted parameters
COEF_UB = [Inf, Inf, 1,Inf];

% Define anonymous function for the fitting
fitfun = @(COEF) (ExpDecay_2Cmp_fun_resTime(COEF, tlist)) - (Y);

% Select Options for the fitting
options = optimset('FunValCheck','off');

% run fitting routibe
[Esp_Coef, resNorm, residuals,exitflag,output,lambda,jacobian] = lsqnonlin(fitfun,COEF0,COEF_LB,COEF_UB,options);

ci = nlparci(Esp_Coef,residuals,'jacobian',jacobian);
Esp_Sigma = (ci(:,2) - ci(:,1))/2;
Esp_Sigma = Esp_Sigma';

if Esp_Coef(2) > Esp_Coef(1)
    tmp = Esp_Coef(1);
    Esp_Coef(1) = Esp_Coef(2);
    Esp_Coef(2) = tmp;
    Esp_Coef(3) = 1 - Esp_Coef(3);
    tmp = Esp_Sigma(1);
    Esp_Sigma(1) = Esp_Sigma(2);
    Esp_Sigma(2) = tmp;
end
    

% Compute output
Esp_Fit(:,1) = data(:,1);
Esp_Fit(:,2) = ExpDecay_2Cmp_fun_resTime(Esp_Coef, data(:,1));

%Compute the indivual components --DB
Esp_Coef_1 = [Esp_Coef(1) Esp_Coef(4)*Esp_Coef(3)]; 
Esp_Fit(:,3) = ExpDecay_fun_resTime(Esp_Coef_1,data(:,1));

Esp_Coef_2 = [Esp_Coef(2) Esp_Coef(4)*(1-Esp_Coef(3))];
Esp_Fit(:,4) = ExpDecay_fun_resTime(Esp_Coef_2,data(:,1));