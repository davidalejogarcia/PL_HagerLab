function out = findParticles(Stack, Threshold, hiBP,windowSz)

nFrames = size(Stack,2);

out = [];

for imageIx = 1:nFrames  
    peaks = pkfnd(Stack(imageIx).data, Threshold, hiBP);
    if ~isempty(peaks)
    centroids = cntrd(Stack(imageIx).data,peaks,windowSz);
    else
    centroids = [0 0 0 0];
    end
    
    if ~isempty(centroids)
        for i = 1: length(centroids(:,1))
        
            x = round(centroids(i,1));
            y = round(centroids(i,2));
            if x ~= 0 && y ~= 0
                centroids(i,5) = Stack(imageIx).data(y,x);
            else
                centroids(i,5) = 0;            
            end
        end
    end
    centroids(:,6) = imageIx;
    out = [out; centroids];
   
     
            
end
    
    
