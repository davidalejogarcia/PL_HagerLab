function [imStack,nFrames] = TIFread(fileName)

imStack = [];

info = imfinfo(fileName);
nFrames = length(info);

for iImg = 1:nFrames
    imStack(iImg).data = imread(fileName,iImg);
    [imStack(iImg).height imStack(iImg).width] = size(imStack(iImg).data);
end



