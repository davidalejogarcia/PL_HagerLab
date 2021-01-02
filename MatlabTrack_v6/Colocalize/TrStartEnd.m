
function OUT = TrStartEnd(Tracks)

nTracks = max(Tracks(:,4));

OUT =zeros(nTracks,6);
for iTrack = 1:nTracks;
    
    idx = find(Tracks(:,4) == iTrack);
    trackTemp = Tracks(idx,:);
    
    AvPositionX = mean(trackTemp(:,1)); 
    AvPositionY = mean(trackTemp(:,2));
    Start = min(trackTemp(:,3));
    End = max(trackTemp(:,3));
    
    OUT(iTrack,:) = [AvPositionX, AvPositionY, ...
        Start, End, End - Start, iTrack];
    
end