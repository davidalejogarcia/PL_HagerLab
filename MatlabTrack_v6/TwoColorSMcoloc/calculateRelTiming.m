function [ArrivalTimes_rel, DepartTimes_rel, coboundTracks1,coboundTracks2] = calculateRelTiming(ImmTracks1,ImmTracks2,ThreshH,Nmin)

%For bound molecules specified in ImmTracks1, and ImmTracks2, find the
%arrival and departure times of molecules in channel 2 relative to channel
%1. ThreshH is the distance that molecules can be separated by and still be
%considered colocalized.

%initialize output variables
ArrivalTimes_rel = [];
DepartTimes_rel = [];
coboundTracks1 = [];
coboundTracks2 = [];

%determine the number of tracks in Channel
nTracks1 = max(ImmTracks1(:,4));
nTracks2 = max(ImmTracks2(:,4));

%calculate the average position of each track
avgTrk1 = zeros(nTracks1,4);
avgTrk2 = zeros(nTracks2,4);

for i = 1:nTracks1
    curT1 = ImmTracks1(ImmTracks1(:,4) == i,:);
    avgTrk1(i,:) = [mean(curT1(:,1)), mean(curT1(:,2)), curT1(1,3), curT1(end,3)];
end
for i = 1:nTracks2
    curT2 = ImmTracks2(ImmTracks2(:,4) == i,:);
    avgTrk2(i,:) = [mean(curT2(:,1)), mean(curT2(:,2)), curT2(1,3), curT2(end,3)];
end

%find colocalized immobile tracks
colocTrks = [];
for i = 1:nTracks1
    curTrk1 = avgTrk1(i,:);
    T1_frames = curTrk1(3):curTrk1(4);
    
    for j = 1:nTracks2
        curTrk2 = avgTrk2(j,:);
        T2_frames = curTrk2(3):curTrk2(4);
        dist = sqrt((curTrk1(1) - curTrk2(1)).^2 + (curTrk1(2) - curTrk2(2)).^2);
        
        overlap_frames = intersect(T1_frames,T2_frames);
        
        if dist <= ThreshH && length(overlap_frames) >= Nmin
            curColoc = [i j];
            colocTrks = [colocTrks; curColoc];
        end
    end
        
    
end

%compile the cobound tracks for channel 1 & 2
for i = 1:size(colocTrks,1)
    trk1 = colocTrks(i,1);
    trk2 = colocTrks(i,2);
    
    
    %compile the relative arrival and departures
    ArrivalTimes_rel = [ArrivalTimes_rel; (avgTrk2(trk2,3) - avgTrk1(trk1,3)), colocTrks(i,:)];
    DepartTimes_rel = [DepartTimes_rel; (avgTrk2(trk2,4) - avgTrk1(trk1,4)), colocTrks(i,:)];
    
end
if ~isempty(colocTrks)

    trk1_all = unique(colocTrks(:,1));
    for i = 1:size(trk1_all,1)

    %compile the cobound tracks for channel 1 & 2
        coboundTracks1 = [coboundTracks1; ImmTracks1(ImmTracks1(:,4) == trk1_all(i),:)];
    end

    trk2_all = unique(colocTrks(:,2));
    for i = 1:size(trk2_all,1)


        coboundTracks2 = [coboundTracks2; ImmTracks2(ImmTracks2(:,4) == trk2_all(i),:)];
    end
end

    
    
    
