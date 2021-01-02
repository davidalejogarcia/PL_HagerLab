function [OUT] = RD_JD_Fun_dz_FreeDb(COEF, PAR)

% Compute the reaction diffusion formula for the Jump distance histogram.
% The formula account for a finite size of the observation slice, according
% to the "absorbing boundary" approximation. The bound popularion is
% allowed to diffuse with a variable diffusion coefficient Db
%
% INPUT
% -------------------------------------------------------------------------
% COEF is a vector with the following structure:
% COEF(1) = D
% COEF(2) = n
% COEF(3) = kon
% COEF(4) = koff
% COEF(5) = DeltaZ
% COEF(6) = Db
%
% PAR is a cell array containing the list of the spatial and temporal
% coordinates for which the Jump histogram has to be calculated.
% PAR{1} = tlist;
% PAR{2} = rlist;
% PAR{e} = sigma;
%
% OUTPUT
% -------------------------------------------------------------------------
% OUT is a matrix with the (time and space dependent) jump distance
% histogram. The structure of OUT is:
%
%        |  OUT(r1, t1)      OUT(r2, t1)    ...      OUT(rmax, t1)   |    
%        |  OUT(r1, t2)      OUT(r2, t2)    ...      OUT(rmax, t2)   |
%  OUT = |      ...              ...        ...           ...        |
%        | OUT(r1, tmax)    OUT(r2, tmax)   ...     OUT(rmax, tmax)  |
%
%--------------------------------------------------------------------------




% Read and prepare input
% ------------------------
tlist = PAR{1};
rlist = PAR{2};
sigma = PAR{3};

D = COEF(1);
n = COEF(2);
kon =  COEF(3);
koff = COEF(4);
deltaZ = COEF(5);
Db = COEF(6);



% Caluculate Jump distance histogram
% -------------------------------------
% The formula for the jump distance histogram is composed by three parts:
% the first part, JDa, represent the molecules that stay bound for the
% whole time t, JDb, the second term JDb is the molecules that stay unbound for the 
% whole time t,the third, JDc that accounts for molecules that associate and 
% dissociate over the inteval t (which need to be integrated numerically)
% Each of the components is represented as a matrix with the following
% structure:
%                |  JDi(r1 , t1)    ....    JDi(rmax, t1)  |
%  JDi(r,t) =    |      ...         ....         ...       |
%                |  JDi(r1, tmax)   ....    JDi(rmax, tmax)|
%


% Calculation of JDa
% JDa is the term accounting for the particles that are bound for the whole
% time t

[rM, tM] = meshgrid(rlist,tlist);

Loc = rM./(2*(sigma^2+Db*tM)).*exp(-rM.^2./(4*(sigma^2+ Db*tM)));      % Spatial dependency
Exp = kon/(kon+koff)*exp(-koff*tM);                                    % Temporal dependency
JDa = Loc.*Exp;                                     

% Calculation of JDb
% JDb is the diffusion solution for the particles that are free for the
% whole time t. It has to be scaled for the probability of staying within
% the observation slice (DeltaZ correction).



% Prepare coefficients for the calculation of the DeltaZ correction
PARdz(:,1) = tlist;
PARdz(1,2) = D;
PARdz(2,2) = 1;

% Call zScale function
zScale = RD_JD_zScaleAB(deltaZ, PARdz, 0);
zScale = zScale*ones(1,length(rlist));



Diff = zScale.*rM./(2*(D*tM + sigma^2)).*...
    exp(-rM.^2./(4*(D*tM + sigma^2))); % diffusive part

Exp = koff/(kon+koff)*exp(-kon*tM); % population decay over time
JDb = Diff.*Exp;

% Calculation of JDc
% JDc can't be calculated matrix-wide, since numerical
% integration is required.

% Set coefficients for numerical integration.
COEFint = [D,kon,koff, sigma, deltaZ, Db];
% initialize Jdc
JDc = zeros(length(tlist), length(rlist));

for i = 1: length(tlist)
    
    IntFun = @(s) RD_JD_Int_dz_BoundDiff(s,tlist(i), rlist, COEFint);
    LowLimit = 0.1*tlist(1);
    HighLimit = 0.99*tlist(i);
    JDc(i,:) = kon*koff / (kon + koff) * quadv(IntFun, LowLimit, HighLimit);


end

%Sum the three components together



OUT = n*(JDa + JDb + JDc)*(rlist(2)-rlist(1));
