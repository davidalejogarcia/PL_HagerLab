function I = gauss_2D_cont_fun(parFit,parMod)

functionFlag = parMod{1};
x = parMod{2};      % Vector with x co-ordinates
y = parMod{3};      % Vector with y co-ordinates

if functionFlag == 1
    I0       = parFit(1);
    sigma_x  = parFit(2);
    sigma_y  = parFit(2);
    bgd      = parFit(3);
    x_center = parMod{4};
    y_center = parMod{5};
elseif functionFlag == 2
    I0       = parFit(1);
    sigma_x  = parFit(2);
    sigma_y  = parFit(2);
    bgd      = parFit(3);
    x_center = parFit(4);
    y_center = parFit(5);  
end

% if sigma_x < 0.25 || sigma_x >10 || sigma_y < 0.25 || sigma_y >10 || x_center < min(x) || x_center > max(x) ||  y_center < min(y) || y_center > max(y) || I0 > 8000
%     I = -1*ones(size(length(y),length(x)));
%     return
% end
    
Ix = exp(-0.5*(x-x_center).^2/sigma_x^2);
Iy = exp(-0.5*(y-y_center).^2/sigma_y^2);

I  = I0*Ix'*Iy + bgd;

