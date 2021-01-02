Wl = 0.700; % Wavelength in microns
NA = 1.45;  % Numerical aperture of the objective
pixelSize = 0.117;
radiusAiry = 2*0.61*Wl/NA/pixelSize;
AiryArea = pi*radiusAiry*radiusAiry;
BGArea = (pi*(radiusAiry+1)*(radiusAiry+1)) - AiryArea;

[process_files,savefold] = uigetfile('*.mat','Select the images to open','','MultiSelect','on');
% savefold = uigetdir(pwd,'Select a save location');

if ~iscell(process_files)
    process_files = {process_files};
end
if process_files{1} ~= 0
    SNR = [];
    int_all = [];
    bg_all = [];
    for i = 1:length(process_files)
        load(fullfile(savefold,process_files{i}));
    %     Particles = Results.Tracking.Tracks(Results.Tracking.Tracks(:,1) > 0,:);
    %     I_bg = true(size(Results.Data.imageStack(1).data));
    %     for j = 1:length(Results.Process.ROIimage)
    %         I_bg(Results.Process.ROIimage{j,:} == 1) = 0;
    %     end
    %     
    %     for j = 1:size(Particles,1)
    %         posX = Particles(j,1);
    %         posY = Particles(j,2);
    %         Time = Particles(j,3);
    %         imCur = Results.Data.imageStack(Time).data;
    %         im_sub = imCur(round(posY-3):round(posY+3),round(posX-3):round(posX+3));
    %         I_peak = max(im_sub(:));
    %         
    %         BG_vals = imCur(I_bg);
    %         int_all = [int_all; double(I_peak)];
    %         bg_all = [bg_all; median(double(BG_vals))];
    %     end

        Tracks = Results.Tracking.Tracks(Results.Tracking.Tracks(:,5) == 1,:);
        tmp = AddAiryIntensity(Tracks,0.104,Results.Data.imageStack);
    %     old_tracks = Results.PreAnalysis.Tracks_um(:,6:7);
    %     Tracks = tmp(:,6:7);
    %     Tracks = Results.PreAnalysis.Tracks_um(:,6:7);
    %     if size(Tracks,1) > 10
    %     Tracks(:,1) = Tracks(:,1)./AiryArea;
    %     Tracks(:,2) = Tracks(:,2)./BGArea;
            int_all = [int_all; double(tmp(:,end-1))];
            bg_all = [bg_all; double(tmp(:,end))];
    %         SNR = [SNR; (int_all - bg_all)./sqrt(int_all)];
    %     end

    end
    SNR = (int_all - bg_all)./sqrt(int_all);
    save([savefold,filesep,'SNR_characteristics.mat'], 'SNR', 'int_all', 'bg_all');
end