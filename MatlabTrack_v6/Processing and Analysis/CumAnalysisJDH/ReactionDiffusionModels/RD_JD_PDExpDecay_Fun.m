function [OUT] = RD_JD_PDExpDecay_Fun(COEF, PAR)

% Compute the pure diffusion formula for the Jump distance histogram
% for a molecule  and diffusing in a 2D space.

%
% INPUT
% -------------------------------------------------------------------------
% COEF is a vector with the following structure:
% COEF(1) = D
% COEF(2) = n
% COEF(3) = sigma
%
% PAR is a cell array containing the list of the spatial and temporal
% coordinates for which the Jump histogram has to be calculated.
% PAR{1} = tlist;
% PAR{2} = rlist;
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

D = COEF(1);
n = COEF(2);
koff = COEF(3);
sigma = COEF(4);

[rM, tM] = meshgrid(rlist,tlist);

% Caluculate Jump distance histogram
% -------------------------------------
% No z-dependency (2D)
zScale = 1;



Diff = zScale.*rM./(2*(D*tM + sigma^2)).*...
    exp(-rM.^2./(4*(D*tM + sigma^2))); % diffusive part

% Calculate Output
OUT = n*exp(-koff*tM).*Diff*(rlist(2)-rlist(1));

