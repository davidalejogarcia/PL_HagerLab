function [filTracks, allFilTracks] = calculateImmobileTracksDim(Tracks, ThresholdL, minFrame, plotFlag)

% calculateImmobileTracks
% --------------------------------------------------
% Find the segments of the tracks that correspond to
% jumps shorter than Threshold.
% ------
% INPUT:
% ------
% tracks: Tracks as produced by the track.m routine. In brief
%         Tracks(:,1) =  x_coordinate (in microns)
%         Tracks(:,2) =  y_coordinate (in microns)
%         Tracks(:,3) =  frame identifier
%         Tracks(:,4) =  particle identifier
% ------
% OUTPUT:
% ------
% filTracks: Immobile track segments:
%       filTracks(:,1) = x-Coordinate (in microns)
%       filTracks(:,2) = y-Coordinate (in microns)
%       filTracks(:,3) = frame identifier
%       filTracks(:,4) = segment identifier
%       filTracks(:,5) = particle identifier
% Different segments with the same particle index belong
% to the same original track


allFilTracks = [];

previous = 0;  % previous = 0 when the previous jump didn't meet
% the filter criterium (jump shorter than 3 * locAcc),
% 1 if it did.
particleIx = 0;

%% Identify particles that Jump Less than ThresholdL from frame to frame;
filTracks = Tracks(Tracks(:,2) >= minFrame,:);
filTracks = filTracks(filTracks(:,3) <= ThreshL,:);
% 
% for i = 1:max(Tracks(:,1));
%     idx = find(Tracks(:,1) == i);
%     
%     trackTemp = Tracks(idx,:);
%     jd_temp = sqrt((trackTemp(2:end,1)- trackTemp(1:end-1,1)).^2 + ...
%         (trackTemp(2:end,2)- trackTemp(1:end-1,2)).^2);
%     
%     for j = 1: length(jd_temp)
%         
%         if jd_temp(j) < ThresholdL && ~previous
%             particleIx = particleIx +1;
%             previous = 1;
%             
%             catTrack = trackTemp ([j j+1],1:3);
%             catTrack(:,4) = particleIx;
%             catTrack(:,[5 6]) = trackTemp ([j j+1],[5 6]);
%             catTrack(:,7) = trackTemp([j j+1], 4);
%             allFilTracks = [allFilTracks ; catTrack];
%             
%             
%         elseif jd_temp(j) < ThresholdL  && previous
%             catTrack = trackTemp (j+1,1:3);
%             catTrack(:,4) = particleIx;
%             catTrack(:,[5 6]) = trackTemp (j+1,[5 6]);
%             catTrack(:,7) = trackTemp (j+1,4);
%             allFilTracks = [allFilTracks ; catTrack];
%             
%         else
%             previous = 0;
%         end
%     end
%     previous = 0;
% end
% 
% %% Discard Tracks shorter than minFrame
% 
% filTracks = discardTracks(allFilTracks,minFrame);
% 
% if  ~isempty(filTracks);
%     %% Identify Tracks that jump less than ThresholdH between minFrame frames
%     
%     
%     for i = 1:max(filTracks(:,4));
%         
%         idx = find(filTracks(:,4)== i);
%         trackTemp = filTracks(idx,:);
%         Test = sqrt((trackTemp(minFrame:end,1)- trackTemp(1:end-minFrame+1,1)).^2 + ...
%             (trackTemp(minFrame:end,2)- trackTemp(1:end-minFrame+1,2)).^2);
%         
%         if max(Test) > ThresholdH;
%             filTracks(idx,:) = [];
%         end
%     end
%     
%     
%     %% Reassign Particle index of filTracks
%     TestAccum = [];
%     Test = filTracks(2:end,4) - filTracks(1:end-1,4) > 0;
%     for i = 1:length(Test);
%         TestAccum(i) = sum(int16(Test(1:i)));
%     end
%     TestAccum = [1, TestAccum + 1];
%     
%     filTracks(:,4) = TestAccum;
%     
% end


if plotFlag
    
    figure;
    hold on
    for i = 1:max(Tracks(:,4));
        idx = find(Tracks(:,4)==i);
        plot(Tracks(idx,1), Tracks(idx,2),'-ok', 'MarkerSize', 5);
    end
    
    for i = 1:max(allFilTracks(:,4));
        idx = find(allFilTracks(:,4)==i);
        plot(allFilTracks(idx,1), allFilTracks(idx,2),'-or', 'MarkerSize', 5);
    end
    
    if  ~isempty(filTracks);
        for i = 1:max(filTracks(:,4));
            idx = find(filTracks(:,4)==i);
            plot(filTracks(idx,1), filTracks(idx,2),'-og', 'MarkerSize', 5);
        end
        
    end
    title('Identified Bound segments');
    xlabel('x [\mum]');
    ylabel('y [\mum]');
end






