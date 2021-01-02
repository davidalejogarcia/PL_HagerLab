function TracksOut = discardTracks(TracksIn,m)
TracksOut = [];
if ~ isempty(TracksIn)
    TracksOut = [];
    ParticleIx = 1;
    
    for i = 1: max(TracksIn(:,4));
        
        idx = find(TracksIn(:,4) == i);
        if length(idx)>=m;
            
            TempTrack = TracksIn(idx,1:end);
            TempTrack(:,4) = ParticleIx;
            TracksOut = [TracksOut; TempTrack];
            ParticleIx = ParticleIx + 1;
        end
        
    end
end