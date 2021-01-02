function TrackLengthHist = calculateTrackLength(tracks, frameTime, minBound)

% This function compute the length of each track and calculate a cumulative
% histogram of the track lengths. The histogram is corrected for
% photobleaching, assuming that an exponential bleaching decay.
%
% Input
% tracks contains the information about the tracks
% in particular track(:, 4) represents the track identifier
% BleachRate is a scalar indicating the photobleaching rate [units of s^-1]
% FrameTime is a scalar indicating the time between two points in the track

if size(tracks,2) == 8
    partInd = 4;
elseif size(tracks,2) == 9
    partInd = 5;
end

nTracks = max(tracks(:,partInd)); 
TrackLength = zeros(nTracks,1);

for i = 1:nTracks;
    idx = find(tracks(:,partInd) == i);  % find the track identified by i
        TrackLength(i) = length(idx);
end


% Compute Number of tracks longer than N-frames
LongestTrack = max(TrackLength);
ShortestTrack = minBound;
TrackLengthHist = zeros(LongestTrack-ShortestTrack + 2,2);
TrackLengthHist(:,1) = ((ShortestTrack-1): LongestTrack)* frameTime;
    
    for  i = 1:nTracks;
        TrackLengthHist(1:TrackLength(i)- ShortestTrack + 1,2) = ...
        TrackLengthHist(1:TrackLength(i) - ShortestTrack + 1,2)+1;
    end;
    
    