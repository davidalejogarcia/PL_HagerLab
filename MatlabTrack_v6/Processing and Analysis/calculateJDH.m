function [tlist, rlist, JDH] = calculateJDH(tracks, nJumps,bin_edges,timeStep, NormalizeFlag, PlotFlag)

% Calculate JDH
% Calculate the Jump histogram distibution from experimental tracks.
% The calculation is calculated for different time lags (DeltaT,
% 2DeltaT,...).
% INPUT:    
% tracks: tracks as produced by the track.m routine. In brief
%         tracks(:,1) =  x_coordinate (in microns)
%         tracks(:,2) =  y_coordinate (in microns)
%         tracks(:,3) =  frame identifier
%         tracks(:,4) =  particle identifier
%
% nJumps: longest time window for jumps calculation 
%         (i.e. nJumps = 10 calculates the jumps for
%          up to 10 frames time interval)
%
% bin_edges: vector containing the edges of the Jump distance histogram
%            bins
%
% timeStep: scalar the time between different frames
%
% Normalize Flag: if 0 each jump is equally contributing to the histogram
%               : if 1 the jump contribute as 1/(LT-n) (with LT being the
%                      length of the track and n being the number of frames
%                      for which the current jump is calculated. 
%                      This allow to prevent overweighting of the Jump 
%                      Histograms at shorter intervals (which produce more
%                      jumps).
%
% Plot Flag: if 0 no plot is produced
%           if 1 a xyz plot is produced
%           if 2 an "heat map" plot is produced
%
% OUTPUT
% tlist = [t1; t2; ...t_maxFrameN)
% rList is a column vector vector with the coordinates corresponding to the centers of
% the histogram bins: rlist = [r1; r2; ...]
% JDH is a matrix with the Jump distance histograms with the following
% structure:
%
%           |  JDH(r1, t1)      JDH(r2, t1)    ...      JDH(rmax, t1)   |    
%           |  JDH(r1, t2)      JDH(r2, t2)    ...      JDH(rmax, t2)   |
%  JDH    = |      ...              ...        ...           ...        |
%           |  JDH(r1, tmax)    JDH(r2, tmax)   ...     JDH(rmax, tmax) |
%
%--------------------------------------------------------------------------


% CALCULATE JUMP HISTOGRAM DISTRIBUTION

% Initialize useful variables
jd = [];                                    % temp variable containing the jumps;
JDH = zeros(nJumps, length(bin_edges));     % Initialize Jump histogram distribution;
nTracks = max(tracks(:,4));                 % number of tracks;
TrackLength = zeros(nTracks, 1);            % Initialize vector containing track length;


for j = 1:nJumps    % for loop on different jump sizes
                    % (1 frame, 2 frames, etc)
    
    for i = 1:nTracks;       % loop on the different tracks.
        
        idx = find(tracks(:,4) == i);       % find the track identified by i
        TrackLength(i) = length(idx);       % compute the length of the i-th track
        
        if  TrackLength(i) >j
            
            jd_temp_x = tracks(idx,1);      % Calculate jumps for the i-th
            jd_temp_y = tracks(idx,2);      % track with the j-th interval
            
            jd_temp = sqrt((jd_temp_x(j+1:end)- jd_temp_x(1:end-j)).^2 + ...
                (jd_temp_y(j+1:end)- jd_temp_y(1:end-j)).^2);
            
            jd =  cat(1,jd,jd_temp);
            
            
            
            if  NormalizeFlag && ~isempty(jd_temp)
                jHist_temp =  histc(jd_temp,  bin_edges);
                % Calculate contribution of the ith track to the histogram at
                % the jth interval
                jHist_temp= jHist_temp/(TrackLength(i)-j); 
                
                % The contribution to the jump histogram is scaled for the
                % number of segments of the track (to avoid overweighting
                % of the jumps at short time)
                

       
                if size(jHist_temp) == [length(bin_edges), 1]
                    jHist_temp = jHist_temp';
                end
                    JDH(j,:) = JDH(j,:) + jHist_temp;
                
                
            end
        end
             
        
    end
    
    if ~NormalizeFlag && ~isempty(jd); % if Normalize flag is set to zero
        
        JDH(j,:) = histc(jd,  bin_edges);           % calculate the JDH for
                                                    % each interval at the
    end                                             % end of the loop on the
                                                    % different tracks
    jd = [];
        
end


JDH = JDH(:,1:end-1);           % Cut the last bin from the histogram
                                % which accounts for the jumps equal
                                % to rmax.
                                
% CALCULATE r AND t;

rlist = bin_edges(1:end-1)+(bin_edges(2)-bin_edges(1))/2;
rlist = rlist';
tlist = (1:nJumps)*timeStep;
tlist = tlist';


% PLOT RESULTS

if PlotFlag == 1
    
    figure;
    hold on;
    
    for t = 1:1:nJumps;
        tplot = tlist(t)*ones(length(rlist), 1);
        plot3(tplot, rlist, JDH(t,:), 'b', 'MarkerSize', 4);
    end
    hold off;
    
    title('Jump Histogram Distribution');
    
    set(gca, 'FontSize', 12);
    view(108,30)
    grid on
    
    xlabel('Time [s]');
    ylabel('Jump Distance [\mum]');
    zlabel('Counts');
        
    
elseif PlotFlag ==2
    
    figure;
    imagesc(rlist, tlist, JDH);
    set(gca, 'FontSize', 12);
    
    xlabel('Jump distance [\mum]');
    ylabel('Time [s]');
    title('Jump Distance Histogram Distribution');
    colorbar('location','EastOutside','XTick',0.5,'XTickLabel','Counts');
    
    axis image
end
        


    









