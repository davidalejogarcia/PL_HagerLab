function I = gauss_step_pix(coord,par)

% Step-wise gaussian function for one pixel
% coord ... coordinates of the pixel
% w     ... FWHM of PSF (instead of the standard deviation sigma)

a = coord(1);
b = coord(2);
c = coord(3);
d = coord(4);

I0 = par{1};
w = par{2};

const = 2*sqrt(log(2))/w;

I = I0*(erf(const*b)-erf(const*a))*(erf(const*d)-erf(const*c));

