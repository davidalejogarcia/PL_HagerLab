function track_length = TrackStartLength

%Displays the scatter plot showing the lenth of tracks vs. when they
%started, and also plots the cumulative distribution of the number of
%tracks started at or after each time-point. The output is a matrix
%containing 4 columns, Column 1: track length, Column 2: image number,
%Column 3: Track number, Column 4: starting frame.


[process_files,pnames] = uigetfile('*.mat','Select the tracked files that you want to analyze','','MultiSelect','on');

if ~iscell(process_files)
    process_files = {process_files};
end
if process_files{1} ~= 0
    Tracks_all = [];
    for i = 1:length(process_files)
        load(fullfile(pnames,process_files{i}));
        Tracks_all = [Tracks_all;Results.Tracking.Tracks];
        Tracks{i} = Results.Tracking.Tracks;
        
    end
    track_length = [];
    for i = 1:length(process_files)
        cur_Track = Tracks{i};
        for j = 1:max(cur_Track(:,4))
            cur = cur_Track(cur_Track(:,4) == j,:);
            cur_len = [size(cur,1),i,j,cur(1,3)];
            track_length = [track_length; cur_len];
        end
    end
    figure; scatter(track_length(:,4),track_length(:,1),'k','filled')
    xlabel('Starting Frame');
    ylabel('Track Length');
    for i = 1:max(track_length(:,4))
        cumhist(i) = length(find(track_length(:,4) > i));
    end
    cumhist_norm = cumhist./cumhist(1);
    
    figure; plot(cumhist)
    xlabel('Frame')
    ylabel Frequency
    title('Tracks starting after Frame i')
end
    
    