function plane = local_background(im_array)

%Calculates the local background around a spot. returns PLANE, which is an
%array the same size as the input image and contains the background to be
%subtracted from each pixel.
% Based on Dan Larson's IDL code.

%David Ball 11/2015

[sizeX, sizeY] = size(im_array);

x_border = double(zeros(2*sizeX,1));
y_border = double(zeros(2*sizeY,1));

x_border(1:sizeX) = im_array(:,1);
x_border(sizeX+1:2*sizeX) = im_array(:,sizeY);

x = double(zeros(2*sizeX,1));
x(1:sizeX) = (1:sizeX)';
x(sizeX+1:2*sizeX) = (1:sizeX)';

y_border(1:sizeY) = im_array(1,:);
yborder(sizeY+1:2*sizeY) = im_array(sizeX,:);

y = double(zeros(2*sizeY,1));
y(1:sizeY) = (1:sizeY)';
y(sizeY+1:2*sizeY) = (1:sizeY)';

delta_x = 2*sizeX*sum(x.^2)-(sum(x)).^2;

a = (1./delta_x)*(sum(x.^2)*sum(x_border)-sum(x)*sum(x.*x_border));
b = (1./delta_x)*(2*sizeX*sum(x.*x_border)-sum(x)*sum(x_border));

delta_y = 2*sizeY*sum(y.^2)-(sum(y)).^2;

c = (1./delta_y)*(sum(y.^2)*sum(y_border)-sum(y)*sum(y.*y_border));
d = (1./delta_y)*(2*sizeY*sum(y.*y_border)-sum(y)*sum(y_border));

offset = (a - d*(sizeY)/2+c - b*(sizeX)/2)/2;


plane = double(zeros(sizeX,sizeY));

for i = 1:sizeX
    for j = 1:sizeY
        plane(i,j) = offset + b*i + d*j;
    end
end
