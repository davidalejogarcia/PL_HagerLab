function [Tracks, Particles] = extendMergeTracks(Tracks,Particles,maxGap)

%Attempts to extend tracks by looking at adjacent frames from each tracks
%extrema, and fitting the last known location to a gaussian. The track will
%only be extended if the fit results in a standard deviation of the
%gaussian less than or equal to sigmaMax. If, while extending the track,
%another track is found that starts at the newly located point, the tracks
%are merged. For best results, the imageStack should be filtered.

NTracks = max(Tracks(:,4));
FrameMax = max(Tracks(:,3));
Trki = 1;
Tracks_tmp = [];






while Trki <= NTracks
    curTrk = Tracks(Tracks(:,4) == Trki,:);
    
    lastFrame = curTrk(end,3);
    TracksNextFr = Tracks(Tracks(:,3) == lastFrame+1,:);
    distFromCur = sqrt((TracksNextFr(:,1) - curTrk(end,1)).^2 + (TracksNextFr(:,2) - curTrk(end,2)).^2);
    blah = 10;
            
        
%     while curTrk(end,3) < FrameMax
%         lkp_X = curTrk(end,1);
%         lkp_Y = curTrk(end,2);
%         %for debugging
%         curImage = imageStack(curTrk(end,3)).data;
%         
%         
%         nextImage = imageStack(curTrk(end,3)+1).data;
%         dim = size(nextImage);
%         dim = [dim(2) dim(1)];
%         % Adjust initial position if the particle is at the edge of the image
%         Pos(1) = min(max(lkp_X, Window + 1), dim(1) - Window);
%         Pos(2) = min(max(lkp_Y, Window + 1), dim(2) - Window);
%         
%         x_sub = (round(Pos(1)) - Window: round(Pos(1)) + Window);
%         y_sub = (round(Pos(2)) - Window: round(Pos(2)) + Window);
%         
%         img_sub = double(nextImage(y_sub,x_sub));
%         % prepare input for the fitting routine
%         coordinates = {x_sub, y_sub};
%         % Fit Gaussian
%         [parFit, ssr] = gauss_fit2D_cov(img_sub,coordinates);
%         
%         if parFit(2) <= sigmaMax
%             TrkPtAdd = [parFit(4), parFit(5), curTrk(end,3) + 1, Trki,0];
%         end
%         
%     end
    
    
    Trki = Trki+1;
end