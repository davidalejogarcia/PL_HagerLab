function results = gmask(im_array,parameters,varargin)

%performs gaussian weighting to extract the total intensity of a
%diffraction-limited spot. 
%Based on Dan Larson's IDL code
%Inputs:
% im_array: image cropped to include the spot and some area surrounding it
% parameters: structure with field psfwidth, used for the size of the
% gaussian
% x0 (optional): center position in X
% y0 (optional): center position in Y
% Specify both x0 & y0 or neither, if they are not specified, gmask will
% attempt to find them
%
%Output:
%results: 1x3 array containing [x0, y0, I], with I being the
%background-subtracted Gaussian-weighted intensity
%
% David Ball 11/2015

[x_dim, y_dim] = size(im_array);
%initialize centroid values to 0
x0 = 0.0;
y0 = 0.0;

F = 1/(sqrt(2)*parameters.psfwidth);

gauss_mask = double(zeros(x_dim,y_dim));

error = 0.0;
results = zeros(1,3);

%perform the local background subtraction
blksubtract = double(im_array) - local_background(im_array);

image = (blksubtract > 0).*blksubtract;

%set the border values to 0
image(1,:) = 0;
image(x_dim,:) = 0;
image(:,1) = 0;
image(:,y_dim) = 0;


[xarr,yarr] = meshgrid(1:x_dim,1:y_dim);

if ~isempty(varargin)
    x0 = varargin{1};
    y0 = varargin{2};
    a = F*(yarr - 0.5 - y0);
    b = F*(yarr + 0.5 - y0);
    c = F*(xarr - 0.5 - x0);
    d = F*(xarr + 0.5 - x0);
    
    gauss_mask = 0.25*(erf(a) - erf(b))*(erf(c) - erf(d));
end

G_total = sum(gauss_mask(:).*gauss_mask(:));
N = sum(image(:).*gauss_mask(:));

photon_number = N/G_total;
results = [x0, y0, photon_number];


