function out = findParticles_3d(Stack, Threshold, hiBP,windowSz)

nFrames = size(Stack,3);

out = [];

for imageIx = 1:nFrames  
    peaks = pkfnd(Stack(:,:,imageIx), Threshold, hiBP);
    if ~isempty(peaks)
    centroids = cntrd(Stack(:,:,imageIx),peaks,windowSz);
    else
    centroids = [0 0 0 0];
    end;
    
    if ~isempty(centroids)
        for i = 1: length(centroids(:,1))
        
            x = round(centroids(i,1));
            y = round(centroids(i,2));
            if x ~= 0 && y ~= 0;
                centroids(i,5) = Stack(y,x,imageIx);
            else
                centroids(i,5) = 0;            
            end
        end
    end
    centroids(:,6) = imageIx;
    out = [out; centroids];
   
     
            
end
    
    
