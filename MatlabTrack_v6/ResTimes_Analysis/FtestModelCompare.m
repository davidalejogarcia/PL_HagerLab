function [h, pval, fstat] = FtestModelCompare(data,model1Fit,model2Fit,p1,p2,varargin)

%Performs a wald test to determine if one model performs significantly
%better. This is a hypothesis test where the null hypothesis is that the
%model with more parameters does not provide a significantly better fit.
% Inputs DATA, MODEL1FIT, and MODEL2FIT should be vectors of the same
% length
% Inputs:
% Data: Experimental data. i.e. values at a given time
% 
% MODEL1FIT: Results of fitting using the first model at the same
% time-points of DATA.
% 
% MODEL2FIT: Results of fitting using the second model at the same
% time-points of DATA.
% 
% P1: Number of parameters in Model 1
%
% P2: Number of parameters in Model 2
% 
% ALPHA (optional): significance level, if not specified the default value
% of 0.05 will be used
%
% Outputs:
% H: Logical result of the hypothesis test. 0 - can't reject the null
% hypothesis at the specified significance level. 1 - can reject the null
% hypothesis
%
% PVAL: p-value of the calculated F-statistic

%Check if alpha is specified
if isempty(varargin)
    alpha = 0.05;
else
    alpha  = varargin{1};
end

%Calculate the SSR for both models
SSR1 = sum((data - model1Fit).^2);
SSR2 = sum((data - model2Fit).^2);

N = length(data);
%calculate the test statistic
if p2 > p1
    if SSR1 > SSR2
        fstat = ((SSR1 - SSR2)/(p2 - p1))/(SSR2/(N-p2));
    else
        fstat = 0;
        
    end
    dof1 = p2 - p1;
    dof2 = N - p2;
elseif p1 > p2
    if SSR2 > SSR1
        fstat = ((SSR2 - SSR1)/(p1 - p2))/(SSR1/(N-p1));
    else
        fstat = 0;
    end
    dof1 = p1 - p2;
    dof2 = N - p1;

    
end
if p1 ~= p2
    pval = 1 - fcdf(fstat,dof1,dof2);
else
    pval = 1;
end

if pval < alpha
    h = 1;
else
    h = 0;
end
