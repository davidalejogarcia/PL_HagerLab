function I_out = registerStack(RegMat)

% [filename1, pathname1] = uigetfile('*.tif','Open image from 488 channel');

[filename2, pathname2] = uigetfile('*.tif','Open image from 647 channel');

% [imageStack488, nImages] = TIFread([pathname1,filesep, filename1]);

[imageStack647,nImages] = TIFread([pathname2,filesep, filename2]);

for i = 1:nImages
    I1 = imageStack647(i).data;
    Roriginal = imref2d(size(I1));
    
    I1_reg = imwarp(I1,RegMat{1},'OutputView',Roriginal);
    
    imwrite(I1_reg,[pathname2,filesep,filename2(1:end-4),'_reg.tif'],'WriteMode','append');
    I_out(:,:,i) = I1_reg;
end