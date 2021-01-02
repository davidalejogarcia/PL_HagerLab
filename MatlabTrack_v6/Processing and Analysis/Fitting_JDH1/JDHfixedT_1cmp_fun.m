function gaussian = JDHfixedT_1cmp_fun(COEF, rlist)

% Calculate the jump distance histogram associated (for a fixed time)
% associated to a single specie diffusing molecules with diffusion
% coefficient D for A jumps.

% Read Input
A=COEF(1);
D = COEF(2);
r = rlist(:,1);
t_bin = rlist(1,2);
dr = r(2) - r(1);

gaussian=dr*r.*A./(2*D*t_bin).*exp(-r.^2/(4*D*t_bin));