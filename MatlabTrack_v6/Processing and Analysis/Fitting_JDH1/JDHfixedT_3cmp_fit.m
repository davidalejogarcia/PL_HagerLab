function [COEF,fit, COEFsigma] = JDHfixedT_3cmp_fit(Hist,jumpTime, D0s)

% JDHfixedT_3cmp_fit
% -------------------------------------------------------------------------
% This function calls the fitting routine for the fit of a jump histogram
% (fixed time). The fitting is performed with the lsqnonlin function from
% from the optimization toolbox. The function is a three-component diffusion
% model
% -------------------------------------------------------------------------
% The output COEF has the resulting parameters of the fit in the order
% [A, D1, D2, f1, SSR];
% The output fit is a two column vector: fit(:,1) =rlist
%                                        fit(:,2) = actual fit



% Read Input
rlist(:,1) = Hist(:,1);
rlist(1,2)=jumpTime;

Y = Hist(:,2);

% define initial values for the parameters
COEF0 = [3*max(Y), D0s(1), D0s(2) D0s(3), 0.33, 0.33];


% define lower boundaries for the fitted parameters
COEF_LB = [0, 0, 0, 0, 0, 0];

% define upper boundaries for the fitted parameters
COEF_UB = [Inf, Inf, Inf, Inf, 1, 1];

% Define anonymous function for the fitting
fitfun = @(COEF) (JDHfixedT_3cmp_fun(COEF, rlist) - Y);



% Define options for the fitting function
options = optimset('FunValCheck','off','Display','off');

% run fitting routibe
[COEF, resNorm, residuals,exitflag,output,lambda,jacobian]= lsqnonlin(fitfun,COEF0,COEF_LB,COEF_UB,options);
ci = nlparci(COEF,residuals,'jacobian',jacobian);
COEFsigma = (ci(:,2) - ci(:,1))/2;
COEFsigma = COEFsigma';

% check if the fit provides acceptable estimates
if min(COEF > 0) == 1 && max(COEF)< Inf
    fit(:,1)=rlist(:,1);
    fit(:,2)=JDHfixedT_3cmp_fun(COEF,rlist);
    
else % Otherwise put the fit to 0
    fit= zeros(length(rlist(:,1)),2);
    disp('THREE COMPONENT FIT of JDH1 FAILED!')
end

SSR = sum((fit(:,2)-Hist(:,2)).^2);

disp('')
disp('--------------------------------------------')
disp('Three component Fit of FJH')
disp(['fitted D1 = ', num2str(COEF(2)), ' mum^2/s'])
disp(['fitted D2 = ', num2str(COEF(3)), ' mum^2/s'])
disp(['fitted D3 = ', num2str(COEF(4)), ' mum^2/s'])
disp(['fitted f1 = ', num2str(COEF(5))])
disp(['fitted f2 = ', num2str(COEF(6))])
disp(['SSR = ', num2str(SSR)]);
disp('--------------------------------------------')
disp('')

COEF(7) = SSR;
