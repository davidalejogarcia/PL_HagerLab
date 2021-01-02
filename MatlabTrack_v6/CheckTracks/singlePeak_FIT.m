function [Fit_Pos,Int,sigma,BG] = singlePeak_FIT(Img, Pos)

% Window size for PSF fitting
dx = 5;
dy = 5;

dim = size(Img);
dim = [dim(2) dim(1)];


% Adjust initial position if the particle is at the edge of the image
Pos(1) = min(max(Pos(1), dx + 1), dim(1) - dx);  
Pos(2) = min(max(Pos(2), dy + 1), dim(2) - dy);

% crop image around the initial position
x_sub = (round(Pos(1)) - dx: round(Pos(1)) + dx);
y_sub = (round(Pos(2)) -dy: round(Pos(2)) + dy);

img_sub = double(Img(y_sub,x_sub));

% prepare input for the fitting routine
 coordinates = {x_sub, y_sub};
 
% Fit Gaussian
 [parFit, ssr] = gauss_fit2D_cov(img_sub,coordinates);
 
% Extract coordinates
 Fit_Pos(1) = parFit(4);
 Fit_Pos(2) = parFit(5);
Int = parFit(1);
sigma = parFit(2);
BG = parFit(3);

