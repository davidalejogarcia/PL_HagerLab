function Integrand = RD_JD_Int_dz_BoundDiff(s,t, rlist, COEF)

% s is a scalar representing the integration variable
% t is a scalar with the maximum time point of integration
% rMax is the maximum jump for which the histogram needs to be calcluated
% rDelta is the space between jump bins
% COEF contains the parameters of kinetic parameters of the JD:
% COEF(:,1) = D;
% COEF(:,2) = kon;
% COEF(:,3) = koff;
% COEF(:,4) = deltaZ;
% -----------------------------------

% Set Parameters
% --------------------------

% Setting the number of elements to be summed together.
% The number of elements to be roughly the product of the total time t by
% the highest between kon and koff (use of the fact that the average number
% of binding events is s*kon and the average number of unbound intervals
% is (t - s)*koff

% Additionally we set a maximum number of elements equal to 1000 and a
% minimum value of 10 
n = max([floor(t*max([COEF(2),COEF(3)])),10]);
n = min([n, 1000]);

n = 1: n;

% Read and prepare input
% ------------------------

D = COEF(1);
kon =  COEF(2);
koff = COEF(3);
sigma = COEF(4);
deltaZ = COEF(5);
Db = COEF(6);



sumLogn = zeros(1, max(n));



% Calculation of the spatial independent term spI(s)
% -----------------------------------------------
% The spatialindependent term is calculated for the first n terms and summed togeter
% before being multipied for the diffusion term


% To avoid calculations with big numbers (producing Inf and NaNs), the
% logarithm of the spatial term is calculated first.

% for the factorial it is useful to calculate the sum of the logarithms of
% n since log(n!) = sum[(x=1 to n)log(x)]

for i = 1:max(n)
    sumLogn(i) = sum(log(1:i));
end;


logAll = log(kon) + log(koff) + log(s) + log(t-s);

spI = (n-1)*(logAll) - 2 * sumLogn + log(n) - ...
(t-s)*koff - s*kon +  + log(s*koff + (t-s)*kon + 2*n);

spI = sum(exp(spI));

% Calculation of the scale factor to account for finite observation thickness
% ----------------------------------------------------------------------------



% Prepare coefficients for the calculation of the DeltaZ correction
PARdz(:,1) = s;
PARdz(1,2) = D;
PARdz(2,2) = 1;

zScale = RD_JD_zScaleAB(deltaZ, PARdz, 1);
zScale = zScale(1);
spI = spI*zScale;


% Calculation of the  spatial dependent term phi(s,r)
% -----------------------------------------------------

% This term is independent of n and is calculated as a vector with the
% following structure:

%               | phi(s, r1)  |      
%               | phi(s, r2)  |       
% phi(s,r) =    |  ...        |
%               | phi(s, rmax)|


phi = rlist/(2* (D * s + Db*(t-s)+ sigma^2))...
    .*exp(-rlist.^2/(4*(D*s + Db*(t-s) + sigma^2)));

Integrand = phi*spI;
