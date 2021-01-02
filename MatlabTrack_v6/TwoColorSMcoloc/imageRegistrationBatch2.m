function imageRegistrationBatch2(inFold,RegMat)

cd(inFold);
sepind = strfind(inFold,filesep);

filenameBase = inFold(sepind(end)+1:end);

% mkdir('aligned');
savefold = [inFold, filesep,'aligned',filesep];
cd('Pos0');

ims = dir('*.tif');
nImages = length(ims);
for i = 1:nImages
    I = imread(ims(i).name);
    if ~isempty(strfind(ims(i).name,'Andor3'))
        frameNum = i/2;
        numstr = num2str(frameNum);
        if frameNum < 10
            numstr = ['0' numstr];
        end
        if frameNum < 100
            numstr = ['0' numstr];
        end
        imwrite(I,[savefold,filenameBase,'_647_t', numstr,'.tif']);
    else
        frameNum = ((i-1)/2) + 1;
        numstr = num2str(frameNum);
        if frameNum < 10
            numstr = ['0' numstr];
        end
        if frameNum < 100
            numstr = ['0' numstr];
        end
        Roriginal =  imref2d(size(I));
        I2_align = imwarp(I,RegMat,'OutputView',Roriginal);
        
        imwrite(I2_align,[savefold,filenameBase,'_488_t', numstr,'.tif']);
    end
end