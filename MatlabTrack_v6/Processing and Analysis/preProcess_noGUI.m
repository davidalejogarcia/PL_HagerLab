function [Tracks_um, NParticles, IntensityHist] = preProcess_noGUI(Tracks,imageStack,...
            Particles, pixelSize, nImages, fileName,ROIpos)


% Import data from the workspace
% handles.Tracks = varargin{1};        % Import tracks.
% handles.Stack = varargin{2};                % Import stacks of images.
% handles.Centroids = varargin{3};     % Import Centroids
% handles.pixelSize = varargin{4};     % Import pixelsize
% handles.NFrames = varargin{5};       % Import number of frames.
% handles.FileName = varargin{6};      %Import FileName.
% handles.ROIpos = varargin{7};
% Find Axes
% handles.TrackAxes = findobj ('Tag', 'Tracks');
% handles.NPartAxes = findobj ('Tag', 'NParticles');
% handles.Hist1Axes = findobj ('Tag', 'Hist1');
% handles.Hist2Axes = findobj ('Tag', 'Hist2');
% 

% Plot All the Tracks.
% axes(handles.TracksAxes);
% SumImage = zeros(length(imageStack(1).data(:,1)),...
%     length(imageStack(1).data(1,:)));
% 
% 
% for imageIx =  1:nImages;
%     SumImage = SumImage +  double(imageStack(imageIx).data);
% end;
% handles.SumImage = SumImage;
% colormap(gray);
% imagesc(SumImage);
% axis image
% hold on;
% hold all;

if size(Tracks,2) < 8
%     for i = 1:1:max(Tracks(:,4))
%         ix = find(Tracks(:,4)==i);
%         plot(Tracks(ix,1),Tracks(ix,2),'LineWidth',1);
%         
%     end;
%     
%     box on;
%     title('Tracks of identified particles');
%     xlabel('X-position [\mum]');
%     ylabel('Y-position [\mum]');
%     
%     hold off;
    if size(Particles,2) > 12
        ROIindx = 13;
    elseif size(Particles,2) > 6 && size(Particles,2) < 13
        ROIindx = 7;
    else
        ROIindx = 0;
    end
    
    if ROIindx > 0
        nROIs = length(ROIpos);
    else
        nROIs = 1;
    end
    
    % Calculate Number of particles for each frame
%     axes(handles.NParticlesAxes);
    NParticles(:,2:nROIs+1) = CalculateNparticles(Particles, nImages,nROIs);
    NParticles(:,1) = 1:length(NParticles(:,2));
   
    
    % Add Airy Intensity to Tracks
    Tracks = AddAiryIntensity(Tracks, ...
        pixelSize, imageStack);
    
    % Calculate histograms of intensities
    IntensityHist = [];
    [IntensityHist(:,2),IntensityHist(:,1)] = hist(Tracks(:,6), 50:50:10000);
    [IntensityHist(:,4),IntensityHist(:,3)] = hist(Tracks(:,7), 50:50:10000);
    [IntensityHist(:,6),IntensityHist(:,5)] = hist(Tracks(:,6)-Tracks(:,7), 50:50:10000);
    
    
    
    
    
    
    
    % Convert the Tracks to microns
    Tracks_um = [];
    Tracks_um(:,1:2) = Tracks(:,1:2)*pixelSize;
    Tracks_um(:,3:7) = Tracks(:,3:7);
    handles.Intensity = Tracks_um(:,6:7);
else
    Intensity = AddAiryIntensityDim(Tracks,pixelSize,imageStack);
    Tracks_um = [];
    Tracks_um(:,1:2) = Tracks(:,1:2);
    Tracks_um(:,3) = Tracks(:,3)*pixelSize;
    Tracks_um(:,4:5) = Tracks(:,4:5);
    Tracks_um(:,6:8) = (Tracks(:,6:8)-1)*pixelSize;
    for i = 1:size(Tracks_um,1)
%         if Tracks_um(i,4) == 1
%             xbox1 = Tracks(i,6);
%             xbox2 = Tracks(i,6) + Tracks(i,3);
%             ybox1 = Tracks(i,7);
%             ybox2 = Tracks(i,8);
%         else
%             ybox1 = Tracks(i,6);
%             ybox2 = Tracks(i,6) + Tracks(i,3);
%             xbox1 = Tracks(i,7);
%             xbox2 = Tracks(i,8);
%         end
%         boxBound = [xbox1, ybox1; xbox1,ybox2; xbox2,ybox2;xbox2,ybox1;xbox1,ybox1];
%         plot(boxBound(:,1),boxBound(:,2));
        
        
    end
%     hold off;
    if size(Particles,2) > 12
        ROIindx = 13;
    elseif size(Particles,2) > 6 && size(Particles,2) < 13
        ROIindx = 7;
    else
        ROIindx = 0;
    end
    
    if ROIindx > 0
        nROIs = max(Particles(:,ROIindx));
    else
        nROIs = 1;
    end
    
    % Calculate Number of particles for each frame
%     axes(handles.NParticlesAxes);
    NParticles(:,2:nROIs+1) = CalculateNparticles(Particles, nImages);
    NParticles(:,1) = 1:length(NParticles(:,2));
%     plot(sum(handles.Nparticles(:,2:end),2), '+k');
%     xlabel('Frame');
%     ylabel('Number of particles');
%     title('Number of detected particles per frame');
%     
    % Calculate histograms of intensities
    IntensityHist = [];
    [IntensityHist(:,2),IntensityHist(:,1)] = hist(Intensity(:,2), 50:50:10000);
    [IntensityHist(:,4),IntensityHist(:,3)] = hist(Intensity(:,3), 50:50:10000);
    [IntensityHist(:,6),IntensityHist(:,5)] = hist(Intensity(:,2)-Intensity(:,3), 50:50:10000);
    
    
    % Plot histograms of Intensity and Background
    
%     axes(handles.Hist1)
%     bar(IntensityHist(:,3),IntensityHist(:,4),'b','EdgeColor', 'none')
%     hold on
%     bar(IntensityHist(:,1),IntensityHist(:,2),'r', 'EdgeColor', 'none')
%     hold off
%     maxPlot = IntensityHist(find(IntensityHist(:,2) ~= 0, 1,'last'),1);
%     xlim([0 maxPlot]);
%     
%     
%     title({'Intensity and Background histograms', 'for tracked particles'});
%     xlabel('Intensity [AU]')
%     ylabel('Counts')
%     legend('Background', 'Particle Intensity');
%     
%     
%     axes(handles.Hist2)
%     bar(IntensityHist(:,5),IntensityHist(:,6),'b')
%     maxPlot = IntensityHist(find(IntensityHist(:,6) ~= 0, 1,'last'),5);
%     xlim([0 maxPlot]);
%     title({'Intensity (background subtracted) histograms', 'for tracked particles'});
%     xlabel('Intensity [AU]')
%     ylabel('Counts')
%     
%     
%     handles.IntensityHist = IntensityHist;
end