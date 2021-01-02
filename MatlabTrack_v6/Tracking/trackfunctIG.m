function [Tracks, TrkPtsAdded, errorcode] = trackfunctIG(Centroid,maxJump,Trackparam);

TrkPtsAdded = {};
if isempty(Centroid)
   errorcode = 5;
   Tracks = [];
return
end
% Tracking
shTr = Trackparam.good;
Trackparam.good = 1;

[Tracks, errorcode] = trackIG(Centroid ,maxJump,Trackparam);

if errorcode ~= 0
    
    return
end
    

%  Deleting 'zeros' lines 
% (correspoding to frames where no track has beenfound

a = find(Tracks(:,1)==0);
Tracks(a,:)= [];


% Fill Gaps (blinking) in tracks
TrkPtsAdded = cell(max(Tracks(:,4)),1);
for i = 1: max(Tracks(:,4))
    pIx = find (Tracks(:,4) == i);
    tempTrack = Tracks(pIx,:);
    diff = tempTrack(2:end,3)-tempTrack(1:end-1,3);
    lIx = find (diff ~=1);
    if  ~isempty(lIx);
        for iGap = 1:length(lIx)
            
            for iFrame = 1:diff(lIx(iGap))-1
                catTrack(:,1) =  tempTrack(lIx(iGap),1) + iFrame*(tempTrack(lIx(iGap)+1,1) - tempTrack(lIx(iGap),1))/diff(lIx(iGap));
                catTrack(:,2) = tempTrack(lIx(iGap),2)+ iFrame*(tempTrack(lIx(iGap)+1,2) - tempTrack(lIx(iGap),2))/diff(lIx(iGap));
                catTrack(:,3) = tempTrack(lIx(iGap),3) + iFrame;%round((tempTrack(lIx+1,3) + tempTrack(lIx,3))/diff(lIx));
                catTrack(:,4) = tempTrack(lIx(iGap),4);         
                Tracks = [Tracks;catTrack];
                %keep a note of which particle points were added
                TrkPtsAdded{i} = [TrkPtsAdded{i}; tempTrack(lIx(iGap),3) + iFrame];
            end
        end
        clear catTrack;
    end
    
%     TrkPtsAdded{i} = (tempTrack(lIx+1,3) + tempTrack(lIx,3))/2;
end



% Re-sort track with the filled gaps
Tracks=sortrows(Tracks, [4 3]);

%remove tracks that are shorter than the specified shortest track
ind = 1;
Tracks2 = [];

for i = 1:max(Tracks(:,4))
    tmp = Tracks(Tracks(:,4) == i,:);
    if size(tmp,1) >= shTr
        tmp(:,4) = ind;
        ind = ind+1;
        Tracks2 = [Tracks2;tmp];
    end
    
end
Tracks = Tracks2;