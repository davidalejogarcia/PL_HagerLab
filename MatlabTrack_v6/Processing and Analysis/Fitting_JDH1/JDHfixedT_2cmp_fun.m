function gaussian = JDHfixedT_2cmp_fun(COEF, rlist)

% Calculate the jump distance histogram associated (for a fixed time)
% associated to two species of diffusing molecules with diffusion
% coefficient D for A jumps.

% Read Input
A=COEF(1);
D1 = COEF(2);
D2 = COEF(3);
f1 = COEF(4);

r = rlist(:,1);
t_bin = rlist(1,2);
dr = r(2) - r(1);

gaussian=dr.*r.*A.*...
    (f1/(2*D1*t_bin).*exp(-r.^2/(4*D1*t_bin)) + ...
 (1 - f1)/(2*D2*t_bin).*exp(-r.^2/(4*D2*t_bin)));