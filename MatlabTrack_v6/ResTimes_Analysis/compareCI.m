function [h, pval, D, CI_tot] = compareCI(val1,val2,varargin)

% statistical test with null-hypothesis that val1 and val2 are the same.
% Inputs val1, and val2 should be 2-element vectors that include the mean 
% and the 95% Confidence interval on the mean. The signicance level can
% also be specified. If not specified, the default 0.05 is used

if isempty(varargin)
    sig = 0.05;
else
    sig = varargin{1};
end
%Extract the mean values
M1 = val1(1);
M2 = val2(1);

%Extract the CI
CI1 = val1(2);
CI2 = val2(2);

%Calculate the difference in the means
D = abs(M1 - M2);

%Combine the CIs
CI_tot = sqrt(CI1.^2 + CI2.^2);

SEM_tot = CI_tot./1.96;

z = D/SEM_tot;
pval = exp(-0.717*z - 0.416*z.^2);

if pval <= sig
    h = 1;
else
    h = 0;
end