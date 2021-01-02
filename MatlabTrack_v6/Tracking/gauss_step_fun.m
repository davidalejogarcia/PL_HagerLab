function I = gauss_step_fun(parFit,parMod)

functionFlag = parMod{1};
x = parMod{2};      % Vector with x co-ordinates
y = parMod{3};      % Vector with y co-ordinates

if functionFlag == 1
    I0       = parFit(1);
    w        = parFit(2);
    bgd      = parFit(3);
    x_center = parMod{4};
    y_center = parMod{5};
elseif functionFlag == 2
    I0       = parFit(1);
    w        = parFit(2);
    bgd      = parFit(3);
    x_center = parFit(4);
    y_center = parFit(5);
elseif functionFlag == 3
    I0       = parFit(1);
    w        = parMod{7};
    bgd      = parMod{6};
    x_center = parFit(2);
    y_center = parFit(3); 
end

if w < 1 || w >10 || x_center < min(x) || x_center > max(x) ||  y_center < min(y) || y_center > max(y) || I0 > 8000
    I = -1*ones(size(length(y),length(x)));
    return
end
    


const = 2*sqrt(log(2))/w;

xx=double(x)-x_center;
yy=double(y)-y_center;

dx = xx(2)-xx(1);
dy = yy(2)-yy(1);

x2 = xx+0.5*dx;
x1 = xx-0.5*dx;

y2 = yy+0.5*dy;
y1 = yy-0.5*dy;

Ix = erf(const*x2)-erf(const*x1);
Iy = erf(const*y2)-erf(const*y1);

I  = I0*Ix'*Iy + bgd;

