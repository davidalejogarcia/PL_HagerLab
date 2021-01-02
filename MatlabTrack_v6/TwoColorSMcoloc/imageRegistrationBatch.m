function imageRegistrationBatch(inFile,RegMat)

%check to see if the output folder already exists
[imageStack, nImages] = TIFread(inFile);
mkdir('aligned');
cd aligned
for i = 1:nImages
    if rem(i,2) == 1
        frameNum = ((i-1)/2) + 1;
        numstr = num2str(frameNum);
        if frameNum < 10
            numstr = ['0' numstr];
        end
        if frameNum < 100
            numstr = ['0' numstr];
        end
        imwrite(imageStack(i).data,[inFile(1:end-20),'647_t', numstr,'.tif']);
    else
        Roriginal =  imref2d(size(imageStack(1).data));
        I2_align = imwarp(imageStack(i).data,RegMat,'OutputView',Roriginal);
        
        frameNum = i/2;
        numstr = num2str(frameNum);
        if frameNum < 10
            numstr = ['0' numstr];
        end
        if frameNum < 100
            numstr = ['0' numstr];
        end
        imwrite(I2_align,[inFile(1:end-20),'488_t', numstr,'.tif']);
    end
end