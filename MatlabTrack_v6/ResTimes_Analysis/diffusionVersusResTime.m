function [DiffCoef,trackLengths] = diffusionVersusResTime(ImmTracks,frameTime)

% calculates the diffusion coefficient for each immobile track, to see if 
% there is any correlation between residence time and diffusion

trackN = unique(ImmTracks(:,4));

DiffCoef = zeros(length(trackN),1);

trackLengths = zeros(length(trackN),1);

for i = 1:length(trackN)
    
    curTrack = ImmTracks(ImmTracks(:,4) == trackN(i),:);
    
    dx = diff(curTrack(:,1));
    dy = diff(curTrack(:,2));
    
    displacement2 = (dx.^2) + (dy.^2);
    d2 = mean(displacement2);
    
    DiffCoef(i) = d2/(4*frameTime);
    trackLengths(i) = frameTime*(size(curTrack,1) - 1);
end
    
