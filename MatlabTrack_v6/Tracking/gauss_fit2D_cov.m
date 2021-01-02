function  [parFit, ssr] = gauss_fit2D_cov(img,coordinates)

x = coordinates{1};
y = coordinates{2};

[px,py] = meshgrid(x,y);
sigma = 1.3;
% Starting values


par0(1) = max(img(:)) -  mean(img(:)); % max intensity
par0(2) = length(x)/8; % sigma
par0(3) = mean(img(:)); % background
par0(4) = x(1) + (x(end) - x(1))/2; % x-center
par0(5) = y(1) + (y(end) - y(1))/2;  % y-center

% boundaries for the fit
LB = [0,0,0,x(1),y(1)];
UB = [65536,10,65536, x(end), y(end)];

% fit options
options = optimset('lsqnonlin');
options.Display = 'none';

% fit
% parFit = lsqnonlin(@(P)objfun(P,px,py,img),par0,LB,UB,options);
parFit = lsqnonlin(@(P)objfun2(P,sigma,px,py,img),par0,LB,UB,options);
residuals = objfun2(parFit,sigma,px,py,img);
ssr = sum(residuals.^2);
parFit(2) = sigma;

% Object function
% --------------------
function residuals = objfun (par, x, y, img)
sigma_matrix = [1/par(2), 0; 0, 1/par(2)];
exponent = [x(:) - par(4), y(:) - par(5)]*sigma_matrix;
model = par(3) + par(1)*exp(-sum(exponent.*exponent,2)/2);
residuals = model - img(:);

 
% Object function
% --------------------
function residuals = objfun2 (par,sigma, x, y, img)
sigma_matrix = [1/sigma, 0; 0, 1/sigma];
exponent = [x(:) - par(4), y(:) - par(5)]*sigma_matrix;
model = par(3) + par(1)*exp(-sum(exponent.*exponent,2)/2);
residuals = model - img(:);





