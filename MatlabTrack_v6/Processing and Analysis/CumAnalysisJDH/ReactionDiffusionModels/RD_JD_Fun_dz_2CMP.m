function [OUT] = RD_JD_Fun_dz_2CMP(COEF, PAR)

% compute the reaction diffusion formula for the Jump distance histogram.
% The formula account for a finite size of the observation slice, according
% to the "Boundary at acquisition time" approximation. It computes the
% reaction diffusion formula with two independent diffusion components
% exchanging with a bound (immobile) component.

% INPUT
% -------------------------------------------------------------------------
% COEF is a vector with the following structure:

% COEF(1) = f1
% COEF(2) = n
% COEF(3) = kon
% COEF(4) = koff
% COEF(6) = deltaZ
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

D1 = PAR{5}(1);
D2 = PAR{5}(2);

f1 = COEF(1);
n = COEF(2);
kon =  COEF(3);
koff = COEF(4);
deltaZ = COEF(5);


COEF_D1 = [D1, f1*n, kon, koff, deltaZ];
COEF_D2 = [D2, (1-f1)*n, kon, koff, deltaZ];

OUT1 =  RD_JD_Fun_dz_BoundDiff(COEF_D1, PAR);
OUT2 =  RD_JD_Fun_dz_BoundDiff(COEF_D2, PAR);

OUT = OUT1 + OUT2;





