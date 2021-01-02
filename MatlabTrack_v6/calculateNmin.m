function [Nmin,prob] = calculateNmin(rmax,T_int,D, varargin)

%Calculates the Nmin to use given the maximum displacement for bound
%particles (rmax), the frame interval (T_int), diffusion coefficient of
%freely diffusing molecules (D), and optionally the threshold probability
%to exclude (P_th). If the threshold probability is not specified, the
%default value of 1% is used.


%Determine if the probability threshold is specified, and set to 1% if it
%isn't

if isempty(varargin)
    P_th = 0.01;
else
    P_th = varargin{1};
end

N = (1:100000)';

expon = (rmax^2)/(4*D*T_int);

P_base = (1 - exp(-expon));

P = P_base.^N;

Nmin = 1 + find(P <= P_th,1,'first');

prob = P(P <= P_th);
prob = prob(1);

