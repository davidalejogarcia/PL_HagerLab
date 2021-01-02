function I_out = merge2colorData(I1,I2,colorID,varargin)

%Combines the input images I1 and I2, into an RGB image, with I1 placed in
%colorID(1) color, and I2 placed in colorID(2) color. Min & max of each image
%can be specified in varargin{3} and varargin{4}. If these are not
%specified, clims default to the min & max of each image.
% colorID should be an 2 element integer vector, with each value specifying
% the color:
% 1-red
% 2-green
% 3-blue
if ~isempty(varargin)
    if length(varargin) > 1
        clims1 = varargin{1};
        clims2 = varargin{2};
    elseif length(varargin) == 1
        clims1 = varargin{1};
        clims2 = [min(I2(:)) max(I2(:))];
    end
else
    clims1 = [min(I1(:)) max(I1(:))];
    clims2 = [min(I2(:)) max(I2(:))];
end

%create a placeholder for the output image
I_out = double(zeros(size(I1,1),size(I1,2),3));

Chan1 = double(I1);
Chan2 = double(I2);

Chan1 = (Chan1 - double(clims1(1)))./(double(clims1(2)) - double(clims1(1)));
Chan2 = (Chan2 - double(clims2(1)))./(double(clims2(2)) - double(clims2(1)));
Chan1(Chan1 < 0) = 0;
Chan2(Chan2 < 0) = 0;

Chan1(Chan1 > 1) = 1;
Chan1(Chan1 > 1) = 1;

I_out(:,:,colorID(1)) = Chan1;
I_out(:,:,colorID(2)) = Chan2;
