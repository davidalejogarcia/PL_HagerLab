function Results = RemoveCheckParticlesNotInCheckTracks(Results)

Tracks = Results.Tracking.Tracks;
CheckTracks = Results.Tracking.CheckTracks;
CheckParticles = Results.Tracking.CheckParticles;
for i = 1:size(Tracks,1)
    Xcoord = Tracks(i,1);
    Ycoord = Tracks(i,2); 
    Tpoint = Tracks(i,3);
    
    idx = find(CheckTracks(:,1) == Xcoord & CheckTracks(:,2) == Ycoord & ...
        CheckTracks(:,3) == Tpoint);
    
    if isempty(idx)
        if size(CheckParticles,2) > 10
            PcoordX = 10;
            PcoordY = 11;
        else
            PcoordX = 1;
            PcoordY = 2;
        end
        PIx = find(CheckParticles(:,PcoordX) == Xcoord & ...
            CheckParticles(:,PcoordY) == Ycoord & ...
            CheckParticles(:,6) == Tpoint);
        if ~isempty(PIx)
            CheckParticles(PIx,:) = [];
        end
        
    end
end
Results.Tracking.CheckParticles = CheckParticles;