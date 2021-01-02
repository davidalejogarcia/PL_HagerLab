function varargout = TwoColorBinding_GUI(varargin)
% TWOCOLORBINDING_GUI MATLAB code for TwoColorBinding_GUI.fig
%      TWOCOLORBINDING_GUI, by itself, creates a new TWOCOLORBINDING_GUI or raises the existing
%      singleton*.
%
%      H = TWOCOLORBINDING_GUI returns the handle to a new TWOCOLORBINDING_GUI or the handle to
%      the existing singleton*.
%
%      TWOCOLORBINDING_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TWOCOLORBINDING_GUI.M with the given input arguments.
%
%      TWOCOLORBINDING_GUI('Property','Value',...) creates a new TWOCOLORBINDING_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TwoColorBinding_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TwoColorBinding_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TwoColorBinding_GUI

% Last Modified by GUIDE v2.5 10-Feb-2017 11:50:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TwoColorBinding_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @TwoColorBinding_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before TwoColorBinding_GUI is made visible.
function TwoColorBinding_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TwoColorBinding_GUI (see VARARGIN)

% Choose default command line output for TwoColorBinding_GUI
handles.output = hObject;
% Import Parameters.
if ~isempty(varargin)
    handles.Flag = varargin{1};
    handles.PathName = varargin{2};
    handles.FileNames = varargin{3};
    AnalysisS = varargin{4};

    handles.Parameters = [];
    handles.Parameters.ThreshL = str2double(AnalysisS{1});
    handles.Parameters.ThreshH = str2double(AnalysisS{2});
    handles.Parameters.minBoundFrames = str2double(AnalysisS{3});
    handles.Parameters.bin = str2double(AnalysisS{4});
    handles.FrameTime = str2double(AnalysisS{5});
    handles.pxSize = str2double(AnalysisS{6});
    handles.Parameters.maxAnalFrame = str2double(AnalysisS{7});
    
    handles.Parameters.Dist2color = str2double(AnalysisS{8});
    handles.Parameters.minInteract = str2double(AnalysisS{9});

    % If only one file has been selected change the fileNames variable
    % to a cell array

    if ~iscell(handles.FileNames)
        handles.FileNames = {handles.FileNames};
    end
    nFiles = length(handles.FileNames);

    % get the tracks of the selected files
    % open each mat-files and retrieve useful information
    ROIChoice = {};
    ROIorClass = -1; % -1: not chosen, 0: ROI name, 1: Class Name
    bgVec = [];
    nROIs = 0;
    
    for i = 1:nFiles

        ROIidx = 0;

        Temp = load([handles.PathName,handles.FileNames{i}],'Results');
        %Below is an attempt to use the image intensities to calculate
        %bleaching
        %     npx(i,:) = Temp.Results.Data.imageStack(1).height.*Temp.Results.Data.imageStack(1).width;
        %
        %     bgVec = [Temp.Results.PreAnalysis.Tracks_um(:,4), Temp.Results.PreAnalysis.Tracks_um(:,7)];
        %     for j = 1:Temp.Results.Data.nImages
        %         BG = bgVec(bgVec(:,1) == j,2);
        %         avgBG = mean(BG);
        %         pxIntesity(i,j) = (mean(Temp.Results.Data.imageStack(j).data(:)) - avgBG)*npx(i);
        %
        %     end
        %Above is an attempt to use the image intensities to calculate
        %bleaching

        if isfield(Temp.Results, 'PreAnalysis1') && ...
                isfield(Temp.Results.PreAnalysis1, 'Tracks_um') && isfield(Temp.Results, 'PreAnalysis2') && ...
                isfield(Temp.Results.PreAnalysis2, 'Tracks_um');
            Tracks1{i} = Temp.Results.PreAnalysis1.Tracks_um;
            if handles.Parameters.maxAnalFrame > 0
                Tracks1{i} = Tracks1{i}(Tracks1{i}(:,3) <= handles.Parameters.maxAnalFrame,:);
            end
            if isfield(Temp.Results.Tracking1,'Particles')
                Particles1{i} = [Temp.Results.Tracking1.Particles(:,10:11), Temp.Results.Tracking1.Particles(:,6)];
            else
                Particles1{i} = [Temp.Results.Tracking1.Centroids(:,1:2), Temp.Results.Tracking1.Particles(:,6)];
            end
            if handles.Parameters.maxAnalFrame > 0
                Particles1{i} = Particles1{i}(Particles1{i}(:,3) <= handles.Parameters.maxAnalFrame,:);
            end
            Particles1{i}(Particles1{i}(:,1) == 0,:) = [];
            Particles1{i}(:,1:2) = (Particles1{i}(:,1:2) - 1).*handles.pxSize;
            if isfield(Temp.Results.Tracking2,'Particles')
                Particles2{i} = [Temp.Results.Tracking2.Particles(:,10:11), Temp.Results.Tracking2.Particles(:,6)];
            else
                Particles2{i} = [Temp.Results.Tracking2.Centroids(:,1:2), Temp.Results.Tracking2.Particles(:,6)];
            end
            if handles.Parameters.maxAnalFrame > 0
                Particles2{i} = Particles2{i}(Particles2{i}(:,3) <= handles.Parameters.maxAnalFrame,:);
            end
            Particles2{i}(Particles2{i}(:,1) == 0,:) = [];
            Particles2{i}(:,1:2) = (Particles2{i}(:,1:2) - 1).*handles.pxSize;
            
            NParticles1{i} = Temp.Results.PreAnalysis1.NParticles;
            Tracks2{i} = Temp.Results.PreAnalysis2.Tracks_um;
            if handles.Parameters.maxAnalFrame > 0
                Tracks2{i} = Tracks2{i}(Tracks2{i}(:,3) <= handles.Parameters.maxAnalFrame,:);
            end
            NParticles2{i} = Temp.Results.PreAnalysis2.NParticles;
            
            if isfield(Temp.Results,'Process1') %check that the strcutures exist
                if isfield(Temp.Results.Process1,'ROIlabel')
                    if isfield(Temp.Results.Process1,'AllROIClasses') && ROIorClass == -1
                        ROIClass =  questdlg('Do you want to separate data based on ROI names or Class Names?','ROI or Class','ROI','Class','Class');
                        if strcmp(ROIClass,'ROI')
                            ROIorClass = 0;
                        else
                            ROIorClass = 1;
                        end
                    elseif ~isfield(Temp.Results.Process1,'AllROIClasses')
                        ROIorClass = 0;
                    end
                    if ROIorClass == 0
                        roiLabels = Temp.Results.Process1.ROIlabel;
                    elseif ROIorClass == 1
                        if ~isfield(Temp.Results.Process1,'AllROIClasses')
                            errordlg(['File does not contain Class data: ,' handles.FileNames{i}]);
                        else
                            roiLabels = Temp.Results.Process1.AllROIClasses;
                        end
                    end
                    if size(roiLabels,1) > 1
                        if ~isempty(ROIChoice) %if a named ROI was selected before check to see if an ROI with that name exists in this dataset
                            if ~strcmp(ROIChoice{1},'All')
                                for j = 1:length(ROIChoice)
                                    for k = 1:length(roiLabels)
                                        if strcmpi(ROIChoice{j},roiLabels{k})
                                            ROIidx = k;
                                            break
                                        end
                                    end
                                end
                            else
                                ROIidx = size(roiLabels,1) + 1; % If all was selected before, just set the index to something greater than the number of ROIs
                            end
                        end

                        if ROIidx == 0 %if nothing has been found automatically, provide the user with a dialog box to select the ROI

                            ROIstring = roiLabels;
                            ROIstring{end+1,1} = 'All';
                            if ROIorClass == 0
                                ROIidx = ROIchooseDlg(ROIstring);
                            else
                                ROIidx = ROIClassChooseDlg(ROIstring);
                            end

                        end
                        %Update the Tracks & NParticles data
                        if ROIidx <= size(roiLabels,1)
                            if ROIorClass == 0
                                if size(Tracks1{i},2) < 8
                                    Tracks1{i} = Tracks1{i}(Tracks1{i}(:,5) == ROIidx,:);
                                    Tracks2{i} = Tracks2{i}(Tracks2{i}(:,5) == ROIidx,:);
                                elseif size(Tracks{i},2) == 8
                                    Tracks1{i} = Tracks1{i}(Tracks1{i}(:,6) == ROIidx,:);
                                    Tracks2{i} = Tracks2{i}(Tracks2{i}(:,6) == ROIidx,:);
                                else
                                    Tracks1{i} = Tracks1{i}(Tracks1{i}(:,9) == ROIidx,:);
                                    Tracks2{i} = Tracks2{i}(Tracks2{i}(:,9) == ROIidx,:);
                                end
                                NParticles1{i} = [NParticles1{i}(:,1) NParticles1{i}(:,ROIidx+1)];
                                NParticles2{i} = [NParticles2{i}(:,1) NParticles2{i}(:,ROIidx+1)];
                                ROIChoice{end+1,1} = Temp.Results.Process.ROIlabel{ROIidx};
                                nROIs = nROIs + 1;
                                
                            else
                                roiAct = [];
                                tracks_tmp1 = Tracks1{i};
                                tracks_tmp2 = Tracks2{i};
                                Tracks1{i} = [];
                                Tracks2{i} = [];
                                for m = 1:size(Temp.Results.Process1.ROIClass)
                                    if strcmpi(roiLabels{ROIidx,:},Temp.Results.Process1.ROIClass{m,:})
                                        if size(Tracks{i},2) < 8
                                            Tracks1{i} = [Tracks1{i}; tracks_tmp1(tracks_tmp1(:,5) == m,:)];
                                            Tracks2{i} = [Tracks2{i}; tracks_tmp2(tracks_tmp2(:,5) == m,:)];
                                        elseif size(Tracks{i},2) == 8
                                            Tracks1{i} = [Tracks1{i}; tracks_tmp1(tracks_tmp1(:,6) == m,:)];
                                            Tracks2{i} = [Tracks2{i}; tracks_tmp2(tracks_tmp2(:,6) == m,:)];
                                        else
                                            Tracks1{i} = [Tracks1{i}; tracks_tmp1(tracks_tmp1(:,9) == m,:)];
                                            Tracks2{i} = [Tracks2{i}; tracks_tmp2(tracks_tmp2(:,9) == m,:)];
                                        end
                                        roiAct = [roiAct;m];
                                    end

                                end
                                tmp1 = NParticles1{i}(:,roiAct+1);
                                tmp2 = NParticles2{i}(:,roiAct+1);
                                nROIs = nROIs + size(tmp1,2);

                                NParticles1{i} = [NParticles1{i}(:,1) sum(tmp1,2)];
                                NParticles2{i} = [NParticles2{i}(:,1) sum(tmp2,2)];
                                ROIChoice{end+1,1} = Temp.Results.Process1.AllROIClasses{ROIidx};
                            end

                        else
                            tmp1 = NParticles1{i}(:,2:size(Temp.Results.Process1.ROIlabel,1)+1);
                            tmp2 = NParticles2{i}(:,2:size(Temp.Results.Process2.ROIlabel,1)+1);
                            nROIs = nROIs + size(tmp1,2);
                            NParticles1{i} = [NParticles1{i}(:,1) sum(tmp1,2)];
                            NParticles2{i} = [NParticles2{i}(:,1) sum(tmp2,2)];
                            ROIChoice{end+1,1} = 'All';
                        end
                    else
                        ROIChoice{1} = roiLabels{1};
                        if ROIorClass == 1
                            tmp1 = NParticles1{i}(:,2:size(NParticles1{i},2));
                            tmp2 = NParticles2{i}(:,2:size(NParticles2{i},2));
                            NParticles1{i} = [NParticles1{i}(:,1) sum(tmp1,2)];
                            NParticles2{i} = [NParticles2{i}(:,1) sum(tmp2,2)];
                            nROIs = nROIs + size(tmp1,2);
                        else
                            nROIs = nROIs + 1;
                        end
                    end
                else
                    ROIChoice{1} = '';
                end
            else
                ROIChoice{1} = '';
            end


            if ~isempty(Tracks1{i}) && ~isempty(Tracks2{i})
                if size(Tracks1{i},2) <= 8
                    [ImmTracks1, Dummy] = calculateImmobileTracks...
                        (Tracks1{i}, handles.Parameters.ThreshL,...
                        handles.Parameters.minBoundFrames, handles.Parameters.ThreshH, 0);
                    [ImmTracks2, Dummy] = calculateImmobileTracks...
                        (Tracks2{i}, handles.Parameters.ThreshL,...
                        handles.Parameters.minBoundFrames, handles.Parameters.ThreshH, 0);
                    [ArrivalTimes_rel{i}, DepartTimes_rel{i}, ...
                        coboundTracks1, coboundTracks2] = calculateRelTiming(...
                        ImmTracks1,ImmTracks2,handles.Parameters.Dist2color,...
                        handles.Parameters.minInteract);
                    ColocPart(i,:) = findColocParticles(Particles1{i},Particles2{i},handles.Parameters.ThreshL);

                    if isempty(coboundTracks1)
                        TrackLengthHist1{i} = [0 0];
                    else
                        
                        coboundTracks1_all{i} = coboundTracks1;
                        % Identify bound molecules for each of the movies;

                        % Calculate survival histogram;
                        TrackLengthHist1{i} = ...
                            calculateTrackLength(coboundTracks1, ...
                            handles.FrameTime,handles.Parameters.minBoundFrames);



                    end
                    if isempty(coboundTracks2)
                        TrackLengthHist2{i} = [0 0];
                    else
                        
                        coboundTracks2_all{i} = coboundTracks2;
                        % Identify bound molecules for each of the movies;

                        % Calculate survival histogram;
                        TrackLengthHist2{i} = ...
                            calculateTrackLength(coboundTracks2, ...
                            handles.FrameTime,handles.Parameters.minBoundFrames);



                    end
                else
                    ImmTracks1 = Tracks1{i}(Tracks1{i}(:,2) >= handles.Parameters.minBoundFrames,:);
                    ImmTracks1 = ImmTracks1(ImmTracks1(:,3) <= handles.Parameters.ThreshL,:);
                    if ~isempty(ImmTracks1)
                        TrackLength1 = ImmTracks1(:,2);
                        LongestTrack1 = max(TrackLength1);
                        ShortestTrack1 = handles.Parameters.minBoundFrames;

                        TrackLengthHist1{i} = zeros(LongestTrack1-ShortestTrack1 + 1,2);
                        TrackLengthHist1{i}(:,1) = (ShortestTrack1-1: LongestTrack1)* handles.FrameTime;
                        %         TrackLengthHist(:,2) = hist(TrackLength,(ShortestTrack:LongestTrack));
                        for  j = 1:length(TrackLength1)
                            TrackLengthHist1{i}(1:TrackLength1(j)- ShortestTrack1 + 1,2) = ...
                                TrackLengthHist1{i}(1:TrackLength1(j) - ShortestTrack1 + 1,2)+1;
                        end
                    end
                    ImmTracks2 = Tracks2{i}(Tracks2{i}(:,2) >= handles.Parameters.minBoundFrames,:);
                    ImmTracks2 = ImmTracks2(ImmTracks2(:,3) <= handles.Parameters.ThreshL,:);
                    if ~isempty(ImmTracks2)
                        TrackLength2 = ImmTracks2(:,2);
                        LongestTrack2 = max(TrackLength2);
                        ShortestTrack2 = handles.Parameters.minBoundFrames;

                        TrackLengthHist2{i} = zeros(LongestTrack2-ShortestTrack2 + 1,2);
                        TrackLengthHist2{i}(:,1) = (ShortestTrack2-1: LongestTrack2)* handles.FrameTime;
                        %         TrackLengthHist(:,2) = hist(TrackLength,(ShortestTrack:LongestTrack));
                        for  j = 1:length(TrackLength2)
                            TrackLengthHist2{i}(1:TrackLength2(j)- ShortestTrack2 + 1,2) = ...
                                TrackLengthHist2{i}(1:TrackLength2(j) - ShortestTrack2 + 1,2)+1;
                        end
                    end
                end
            end


        else
            errordlg(['Some of the files have not been preprocessed: ',handles.PathName,handles.FileNames{i}]);
            return
        end
    end
    handles.coboundTracks1_all = coboundTracks1_all;
    handles.coboundTracks2_all = coboundTracks2_all;
    handles.nROIs = nROIs;
    %Calculate the avg image intensity over time
    % avgInt = sum(pxIntesity)./sum(npx);

    % Calculate the total number of particles;

    % find the data containing the lowest number of frames;
%     TimePoints = [];
    ArrivalTimes_rel_all = [];
    DepartTimes_rel_all = [];
    TrkNum_ArrDep = [];
    for i = 1:nFiles
        if ~isempty(ArrivalTimes_rel{i})
            ArrivalTimes_rel_all = [ArrivalTimes_rel_all; ArrivalTimes_rel{i}(:,1)];
            DepartTimes_rel_all = [DepartTimes_rel_all; DepartTimes_rel{i}(:,1)];
            CellNum = i*ones(size(ArrivalTimes_rel{i},1),1);
            TrkNum_ArrDep = [TrkNum_ArrDep; CellNum, ArrivalTimes_rel{i}(:,2:3)];
        end
        
    end
    handles.ArrivalTimes_rel = ArrivalTimes_rel_all;
    handles.DepartTimes_rel = DepartTimes_rel_all;
    handles.TrkNum_ArrDep = TrkNum_ArrDep;
    
    
    if isempty(ArrivalTimes_rel_all) && isempty(DepartTimes_rel_all)
        msgbox('No Tracks are colocalized based on the specified parameters','No Colocalization');
        close('TwoColorBinding_GUI');
        return;
    end
    ArrivalTimes_rel_all = ArrivalTimes_rel_all*handles.FrameTime;
    DepartTimes_rel_all = DepartTimes_rel_all*handles.FrameTime;
    
    Arrive_mean = mean(ArrivalTimes_rel_all);
    Arrive_median = median(ArrivalTimes_rel_all);
    Arrive_std = std(ArrivalTimes_rel_all);
    Arrive_mode = mode(ArrivalTimes_rel_all);
    
    Depart_mean = mean(DepartTimes_rel_all);
    Depart_median = median(DepartTimes_rel_all);
    Depart_std = std(DepartTimes_rel_all);
    Depart_mode = mode(DepartTimes_rel_all);
    
    bin_min = min(min(ArrivalTimes_rel_all),min(DepartTimes_rel_all));
    bin_max = max(max(ArrivalTimes_rel_all),max(DepartTimes_rel_all));
    
    MaxMax = max(bin_max,abs(bin_min));
    
    bins = -1*MaxMax:handles.FrameTime:MaxMax;
    if length(bins) == 1
        bins = 1;
    end
    [handles.ArrivalHist,x] = hist(ArrivalTimes_rel_all,bins);
    handles.ArrivalHist = [x', handles.ArrivalHist'];
    [handles.DepartHist,x] = hist(DepartTimes_rel_all,bins);
    handles.DepartHist = [x', handles.DepartHist'];
    
    axes(handles.ArrivalAxes);
    bar(handles.ArrivalHist(:,1),handles.ArrivalHist(:,2));
    title({['Relative Arrival - mean = ', num2str(Arrive_mean,2), ', std. dev. = ', num2str(Arrive_std,2)],...
        ['median = ', num2str(Arrive_median,2), ', mode = ',num2str(Arrive_mode,2), ', N = ',num2str(size(ArrivalTimes_rel_all,1))]});
    xlabel('Time (s)');
    ylabel('Frequency');
    
    
    axes(handles.DepartAxes);
    bar(handles.DepartHist(:,1),handles.DepartHist(:,2));
    title({['Relative Departure - mean = ', num2str(Depart_mean,2), ', std. dev. = ', num2str(Depart_std,2)],...
        ['median = ', num2str(Depart_median,2), ', mode = ',num2str(Depart_mode,2), ', N = ',num2str(size(DepartTimes_rel_all,1))]});
    xlabel('Time (s)');
    ylabel('Frequency');
    xlim = get(gca,'Xlim');
    
    %Generate plot of Residence times versus arrival & departures
    CoboundTrackLengthsArrDep = zeros(size(TrkNum_ArrDep,1),2);
    for i = 1:size(TrkNum_ArrDep,1)
        curCell1 = handles.coboundTracks1_all{TrkNum_ArrDep(i,1)};
        curTrk1 = curCell1(curCell1(:,4) == TrkNum_ArrDep(i,2),:);
        
        curCell2 = handles.coboundTracks2_all{TrkNum_ArrDep(i,1)};
        curTrk2 = curCell2(curCell2(:,4) == TrkNum_ArrDep(i,3),:);
        
        CoboundTrackLengthsArrDep(i,1) = (size(curTrk1,1)-1)*handles.FrameTime;
        CoboundTrackLengthsArrDep(i,2) = (size(curTrk2,1)-1)*handles.FrameTime;
        
    end
    axes(handles.ResVsArrAxes)
    scatter(ArrivalTimes_rel_all,CoboundTrackLengthsArrDep(:,1),'r');
    hold on
    scatter(ArrivalTimes_rel_all,CoboundTrackLengthsArrDep(:,2),'g');
    hold off
    set(gca,'Xlim',xlim);
    xlabel('Relative Arrival of Channel 2 (s)')
    ylabel('Residence Time (s)');
    legend('Channel 1','Channel 2');
    
    
    axes(handles.ResVsDepAxes)
    scatter(DepartTimes_rel_all,CoboundTrackLengthsArrDep(:,1),'r');
    hold on
    scatter(DepartTimes_rel_all,CoboundTrackLengthsArrDep(:,2),'g');
    hold off
    set(gca,'Xlim',xlim);
    xlabel('Relative Departure of Channel 2 (s)')
    legend('Channel 1','Channel 2');
    
    nMax1 = length(NParticles1{1}(:,1));
    TimePoints1 = NParticles1{1}(:,1);
    for i = 1:nFiles
        
        if size(NParticles1{i},1) < nMax1

            TimePoints1 = NParticles1{i}(:,1);
            nMax1 = length(TimePoints1);
        end
    end
    
    nMax2 = length(NParticles2{1}(:,1));
    TimePoints2 = NParticles2{1}(:,1);
    for i = 1:nFiles
        if isempty(NParticles2{i})
            blag = 0;
        end
        if size(NParticles2{i},1) < nMax2

            TimePoints2 = NParticles2{i}(:,1);
            nMax2 = length(TimePoints2);
        end
    end

    if isempty(TimePoints1) || isempty(TimePoints2)
        errordlg('Too few bound particles to produce an histogram');
        close(TwoColorBinding_GUI);
        return
        
    end

    % Accumulate the histogram
    Hist_Matrix1 = zeros(nMax1,nFiles);
    Hist_Matrix2 = zeros(nMax2,nFiles);

    for i = 1:nFiles
        Hist_Matrix1(1:nMax1,i) = NParticles1{i}(1:nMax1,2);
        Hist_Matrix2(1:nMax2,i) = NParticles2{i}(1:nMax2,2);
        
    end

    CumNParticles1(:,1) = (TimePoints1-1) * handles.FrameTime;
    CumNParticles1(:,2) = sum(Hist_Matrix1, 2);
    CumNParticles2(:,1) = (TimePoints2-1) * handles.FrameTime;
    CumNParticles2(:,2) = sum(Hist_Matrix2, 2);

    %Calculate the avg intensity


    % Accumulate the different survival histograms;

    % find the histogram containing the longest track
    nMax1 = 1;
    TimePoints1 = [];
    totalTracks1 = 0;
    nMax2 = 1;
    TimePoints2 = [];
    totalTracks2 = 0;
    for i = 1:nFiles
        trInd1 = unique(Tracks1{i}(:,4));
        trInd2 = unique(Tracks2{i}(:,4));
        totalTracks1 = totalTracks1 + length(trInd1);
        totalTracks2 = totalTracks2 + length(trInd2);
        if isempty(TrackLengthHist1{i})
            blag = 0;
        end
        if size(TrackLengthHist1{i},1) > nMax1;

            TimePoints1 = TrackLengthHist1{i}(:,1);
            nMax1 = length(TimePoints1);
        end
        if size(TrackLengthHist2{i},1) > nMax2;

            TimePoints2 = TrackLengthHist2{i}(:,1);
            nMax2 = length(TimePoints2);
        end
    end
    handles.totalTracks1 = totalTracks1;
    handles.totalTracks2 = totalTracks2;
    if isempty(TimePoints1) || isempty(TimePoints2)
        errordlg('Too few bound particles to produce an histogram');
        return
    end


    % Accumulate the histograms
    Hist_Matrix1 = zeros(nMax1,nFiles);
    Hist_Matrix2 = zeros(nMax2,nFiles);

    for i = 1:nFiles
        if ~isempty(TrackLengthHist1{i})
            lengthHist1 = length(TrackLengthHist1{i}(:,1));
            Hist_Matrix1(1:lengthHist1,i) = TrackLengthHist1{i}(:,2);
        end
        if ~isempty(TrackLengthHist2{i})
            lengthHist2 = length(TrackLengthHist2{i}(:,1));
            Hist_Matrix2(1:lengthHist2,i) = TrackLengthHist2{i}(:,2);
        end
    end

    CumTrackLengthHist1(:,1) = TimePoints1;
    CumTrackLengthHist1(:,2) = sum(Hist_Matrix1, 2);
    CumTrackLengthHist2(:,1) = TimePoints2;
    CumTrackLengthHist2(:,2) = sum(Hist_Matrix2, 2);
    % Copy Histogram to handles
    handles.CumHist1 = CumTrackLengthHist1;
    handles.CumHist2 = CumTrackLengthHist2;


    % Calculate the fraction of bound molecules.

    % First calculate the residence time histogram
    ResTimeHist1(:,1)= CumTrackLengthHist1(:,1);
    ResTimeHist1(:,2) = ...
        [(CumTrackLengthHist1(1:end-1,2) - CumTrackLengthHist1(2:end,2));CumTrackLengthHist1(end,2)];
    ResTimeHist2(:,1)= CumTrackLengthHist2(:,1);
    ResTimeHist2(:,2) = ...
        [(CumTrackLengthHist2(1:end-1,2) - CumTrackLengthHist2(2:end,2));CumTrackLengthHist2(end,2)];
    
    % Then calculate the total number of bound spots
    TotalBoundMolecules1 = sum(ResTimeHist1(:,2).*ResTimeHist1(:,1)/handles.FrameTime);
    handles.BoundMolecules1 = TotalBoundMolecules1;
    
    TotalBoundMolecules2 = sum(ResTimeHist2(:,2).*ResTimeHist2(:,1)/handles.FrameTime);
    handles.BoundMolecules2 = TotalBoundMolecules2;
    % and the total number of spots
    % TotalMolecules = sum(CumNParticles(:,2));
    TotalMolecules1 = 0;
    TotalMolecules2 = 0;
    for i = 1:nFiles
        TotalMolecules1 = TotalMolecules1 + sum(NParticles1{i}(:,2));
        TotalMolecules2 = TotalMolecules2 + sum(NParticles2{i}(:,2));
    end
    TotalColocMolecules = sum(ColocPart);
    
    handles.TotalMolecules1 = TotalMolecules1;
    handles.TotalMolecules2 = TotalMolecules2;
    % Finally divide the two.
%     PartialBoundFraction1 = TotalBoundMolecules1/TotalMolecules1;
    PartialBoundFraction1 = TotalBoundMolecules1/TotalColocMolecules;
    BFerror1 = (sqrt(TotalBoundMolecules1)/TotalBoundMolecules1 + ...
        sqrt(TotalColocMolecules)/TotalColocMolecules)*PartialBoundFraction1;
    
%     PartialBoundFraction2 = TotalBoundMolecules2/TotalMolecules2;
    PartialBoundFraction2 = TotalBoundMolecules2/TotalColocMolecules;
    BFerror2 = (sqrt(TotalBoundMolecules2)/TotalBoundMolecules2 + ...
        sqrt(TotalColocMolecules)/TotalColocMolecules)*PartialBoundFraction2;

    % NOTE: THIS IS ONLY THE PARTIAL BOUND FRACTION
    % BECAUSE WE ARE LOOKING ONLY AT PARTICLES BOUND
    % FOR MORE THAN Nmin. We HAVE TO FIT THE HISTOGRAM OF RESIDENCE
    % TIMES AND EXTRAPOLATE THE TRUE VALUE OF THE BOUND FRACTION

    % Photobleaching correction




    % Fit a double exponential to the photobleaching curve
    disp(' ')
    disp('____________________________________')
    disp('Estimating bleaching characteristics, Channel 1')
    disp('____________________________________')

    [BleachRates1, Dummy, CumNParticles1(:,3:6)] =...
        ExpDecay_2Cmp_fit(CumNParticles1, [1 0.1]);
    disp(['Bleach Rate 1: ', num2str(BleachRates1(1), 3), ' s^-1'])
    disp(['Bleach Rate 2: ', num2str(BleachRates1(2), 3), ' s^-1'])
    disp(['Fraction 1: ', num2str(BleachRates1(3), 3)])
    disp('____________________________________')
    disp(' ')

    CumNParticles1 = CumNParticles1(:,1:4);
    % Plot the photobleaching curve
%     PhotobAxes = findobj ('Tag', 'photob_axes');
%     axes(PhotobAxes);
% 
% 
% 
% 
% 
% 
%     plot(CumNParticles(:,1), CumNParticles(:,2)/max(CumNParticles(:,4)), '+', 'MarkerEdgeColor', [0.5 0.5 0.5]');
%     hold on
%     plot(CumTrackLengthHist(:,1), CumTrackLengthHist(:,2)/max(CumTrackLengthHist(:,2)),'ok');
%     plot(CumNParticles(:,3), CumNParticles(:,4)/max(CumNParticles(:,4)),'r');
%     hold off;
% 
%     box on;
% 
%     title({'Photobleaching kinetics:', ['k_{b1} = ', num2str(BleachRates(1), 3),...
%         ' s^{-1}, k_{b2} = ', num2str(BleachRates(2), 3), ' s^{-1}, f_1 = ',  ...
%         num2str(BleachRates(3), 3)]})
% 
%     hold off;
%     legend('Photobleaching Decay', 'Bound molecules Decay');
% 
%     xlabel('Time [s]');
%     ylabel('Normalized Counts');
%     % Copy CumNparticles to handles
    handles.NParticles1 = CumNParticles1;
    handles.BleachRates1 = BleachRates1;
    handles.ColocParticles = TotalColocMolecules;
    
    % Fit a double exponential to the photobleaching curve
    disp(' ')
    disp('____________________________________')
    disp('Estimating bleaching characteristics, Channel 2')
    disp('____________________________________')

    [BleachRates2, Dummy, CumNParticles2(:,3:6)] =...
        ExpDecay_2Cmp_fit(CumNParticles2, [1 0.1]);
    disp(['Bleach Rate 1: ', num2str(BleachRates2(1), 3), ' s^-1'])
    disp(['Bleach Rate 2: ', num2str(BleachRates2(2), 3), ' s^-1'])
    disp(['Fraction 1: ', num2str(BleachRates2(3), 3)])
    disp('____________________________________')
    disp(' ')

    CumNParticles2 = CumNParticles2(:,1:4);
    handles.NParticles2 = CumNParticles2;
    handles.BleachRates2 = BleachRates2;

    CumTrackLengthHist_PB1 = CumTrackLengthHist1;
    CumTrackLengthHist_PB1(:,2) = CumTrackLengthHist_PB1(:,2)./ ...
        (BleachRates1(3)*exp(-BleachRates1(1).* CumTrackLengthHist_PB1(:,1)) + ...
        (1-BleachRates1(3))*exp(-BleachRates1(2).* CumTrackLengthHist_PB1(:,1)));
    ResTimeHist_PB1(:,1)= CumTrackLengthHist_PB1(:,1);
    ResTimeHist_PB1(:,2) = ...
        [(CumTrackLengthHist_PB1(1:end-1,2) - CumTrackLengthHist_PB1(2:end,2));CumTrackLengthHist_PB1(end,2)];
    
    CumTrackLengthHist_PB2 = CumTrackLengthHist2;
    CumTrackLengthHist_PB2(:,2) = CumTrackLengthHist_PB2(:,2)./ ...
        (BleachRates2(3)*exp(-BleachRates2(1).* CumTrackLengthHist_PB2(:,1)) + ...
        (1-BleachRates2(3))*exp(-BleachRates2(2).* CumTrackLengthHist_PB2(:,1)));
    ResTimeHist_PB2(:,1)= CumTrackLengthHist_PB2(:,1);
    ResTimeHist_PB2(:,2) = ...
        [(CumTrackLengthHist_PB2(1:end-1,2) - CumTrackLengthHist_PB2(2:end,2));CumTrackLengthHist_PB2(end,2)];

    %only use values where there are particles with that residence time for curve-fitting
    ind2fit1 = find(ResTimeHist1(:,2) > 0);
    ind2fit2 = find(ResTimeHist2(:,2) > 0);


    % Ask if you want to correct the histogram for photobleaching;
    answer = questdlg('Do you want to correct the histogram for photobleaching?');

    switch answer
        case 'Yes'
            % Recalculate the survival probability corrected for photobleaching
    %         CumTrackLengthHist(:,2) = CumTrackLengthHist(:,2)./ ...
    %             (BleachRates(3)*exp(-BleachRates(1).* CumTrackLengthHist(:,1)) + ...
    %             (1-BleachRates(3))*exp(-BleachRates(2).* CumTrackLengthHist(:,1)));
    %         % Recalculate the residence time distribution corrected for
    %         % photobleaching
    %         ResTimeHist(:,1)= CumTrackLengthHist(:,1);
    %         ResTimeHist(:,2) = ...
    %             [(CumTrackLengthHist(1:end-1,2) - CumTrackLengthHist(2:end,2));CumTrackLengthHist(end,2)];
            set(handles.PBCorrect,'Value',1);
        case 'No'
            set(handles.PBCorrect,'Value',0);

    end
    deltaT1 = ResTimeHist1(2,1) - ResTimeHist1(1,1);
    deltaT2 = ResTimeHist2(2,1) - ResTimeHist2(1,1);
    % Normalize the Survival probability for the partial bound fraction
    % if size(Tracks{1},2) < 8
    CumTrackLengthHist1(:,2) = CumTrackLengthHist1(:,2)/CumTrackLengthHist1(1,2)...
        .*PartialBoundFraction1;
    CumTrackLengthHist2(:,2) = CumTrackLengthHist2(:,2)/CumTrackLengthHist2(1,2)...
        .*PartialBoundFraction2;

    CumTrackLengthHist_PB1(:,2) = CumTrackLengthHist_PB1(:,2)/CumTrackLengthHist_PB1(1,2)...
        .*PartialBoundFraction1;
    CumTrackLengthHist_PB2(:,2) = CumTrackLengthHist_PB2(:,2)/CumTrackLengthHist_PB2(1,2)...
        .*PartialBoundFraction2;
    
    % Normalize the residence time histogram for the partial bound fraction
    ResTimeHist1(:,2) = ResTimeHist1(:,2)/(deltaT1*sum(ResTimeHist1(:,2)))*...
        PartialBoundFraction1;
    ResTimeHist2(:,2) = ResTimeHist2(:,2)/(deltaT2*sum(ResTimeHist2(:,2)))*...
        PartialBoundFraction2;

    ResTimeHist_PB1(:,2) = ResTimeHist_PB1(:,2)/(deltaT1*sum(ResTimeHist_PB1(:,2)))*...
        PartialBoundFraction1;
    ResTimeHist_PB2(:,2) = ResTimeHist_PB2(:,2)/(deltaT2*sum(ResTimeHist_PB2(:,2)))*...
        PartialBoundFraction2;
    % end


    % handles.Hist.Surv(:,1:2) = CumTrackLengthHist;
    % 
    % handles.Hist.Surv_PB(:,1:2) = CumTrackLengthHist_PB;
    % 
    % handles.Hist.Res(:,1:2) = ResTimeHist;
    % 
    % handles.Hist.Res_PB(:,1:2) = ResTimeHist_PB;
    handles.PartBoundFrac1 = PartialBoundFraction1;
    handles.PartBoundFrac2 = PartialBoundFraction2;
    handles.BFerror1 = BFerror1;
    handles.BFerror2 = BFerror2;

    % If the user wanted to analyze the residence time, convert the cumulative
    % histogram to a residence time histogram.

    if handles.Flag == 1;
        set(handles.What2Plot,'Value',1);
    else
        set(handles.What2Plot,'Value',2);
    end

        % Bin the histogram according to the selected settings;
        nBins = handles.Parameters.bin;
        if nBins ~= 0
            binFactor1 = floor(length(ResTimeHist1(:,1))/nBins);
            binFactor2 = floor(length(ResTimeHist2(:,1))/nBins);
    %     binFactor = handles.Parameters.bin + 1;
    %     nBins = floor(length(ResTimeHist(:,1))/binFactor);
    %     if strcmp(answer,'Yes')
            TimeBins1 = reshape(ResTimeHist1(1:nBins*binFactor1,1), binFactor1, nBins);
            TimeBins2 = reshape(ResTimeHist2(1:nBins*binFactor2,1), binFactor2, nBins);
            CountsBins1 = reshape(ResTimeHist1(1:nBins*binFactor1,2), binFactor1, nBins);
            CountsBins2 = reshape(ResTimeHist2(1:nBins*binFactor2,2), binFactor2, nBins);
    %     else
            TimeBins_PB1 = reshape(ResTimeHist_PB1(1:nBins*binFactor1,1), binFactor1, nBins);
            TimeBins_PB2 = reshape(ResTimeHist_PB2(1:nBins*binFactor2,1), binFactor2, nBins);
            CountsBins_PB1 = reshape(ResTimeHist_PB1(1:nBins*binFactor1,2), binFactor1, nBins);
            CountsBins_PB2 = reshape(ResTimeHist_PB2(1:nBins*binFactor2,2), binFactor2, nBins);
    %     end


            resTimeHist_Binned1 = [];
            resTimeHist_Binned1(:,1) = mean(TimeBins1);
            resTimeHist_Binned1(:,2) = sum(CountsBins1);
            resTimeHist_Binned2 = [];
            resTimeHist_Binned2(:,1) = mean(TimeBins2);
            resTimeHist_Binned2(:,2) = sum(CountsBins2);
            
            resTimeHist_Binned_PB1 = [];
            resTimeHist_Binned_PB1(:,1) = mean(TimeBins_PB1);
            resTimeHist_Binned_PB1(:,2) = sum(CountsBins_PB1);
            resTimeHist_Binned_PB2 = [];
            resTimeHist_Binned_PB2(:,1) = mean(TimeBins_PB2);
            resTimeHist_Binned_PB2(:,2) = sum(CountsBins_PB2);
        else
            resTimeHist_Binned1 = ResTimeHist1;
            resTimeHist_Binned2 = ResTimeHist2;
            resTimeHist_Binned_PB1 = ResTimeHist_PB1;
            resTimeHist_Binned_PB2 = ResTimeHist_PB2;
        end

        % Remove values < 0 from the histogram for fitting only
        idx1 = find(resTimeHist_Binned1(:,2) < 0);
        resTimeHist_Binned1(idx1,2) = 0;
        idx1 = find(resTimeHist_Binned1(:,2) <= 0);
        resTimeHist_Binned1_2 = resTimeHist_Binned1;
        resTimeHist_Binned1_2(idx1,:) = [];
        % Fit the binned residence time histogram with exponentials;
        [fitpar1,espSigma1, ~]= ExpDecay_fit_resTime(resTimeHist_Binned1_2, 1 );
        fit1 = [resTimeHist_Binned1(:,1),ExpDecay_fun_resTime(fitpar1,resTimeHist_Binned1(:,1))];
        [fitpar1_2,espSigma1_2, ~]= ExpDecay_2Cmp_fit_resTime(resTimeHist_Binned1_2, [0.1 0.01]);

        Esp_Coef1_1 = [fitpar1_2(1) fitpar1_2(4)*fitpar1_2(3)];
        Esp_Fit1(:,1) = ExpDecay_fun_resTime(Esp_Coef1_1,resTimeHist_Binned1(:,1));

        Esp_Coef1_2 = [fitpar1_2(2) fitpar1_2(4)*(1 - fitpar1_2(3))];
        Esp_Fit1(:,2) = ExpDecay_fun_resTime(Esp_Coef1_2,resTimeHist_Binned1(:,1));

        fit2 = [resTimeHist_Binned1(:,1),ExpDecay_2Cmp_fun_resTime(fitpar1_2,resTimeHist_Binned1(:,1)),Esp_Fit1];

        handles.Hist.Res1(:, 1:2) = resTimeHist_Binned1;
        handles.Hist.Res1(:,3) = fit1(:,2);
        handles.Hist.Res1(:,4) = fit2(:,2);
        handles.Hist.Res1(:,5) = fit2(:,3);
        handles.Hist.Res1(:,6) = fit2(:,4);

        [~, pvalRes1,FstatRes1] = FtestModelCompare(resTimeHist_Binned1(:,2),fit1(:,2),fit2(:,2),2,4);

        handles.FitPar.Res1 = [fitpar1, fitpar1_2, PartialBoundFraction1, pvalRes1; espSigma1, espSigma1_2, BFerror1,FstatRes1];
        idx = find(resTimeHist_Binned_PB1(:,2) < 0);
        resTimeHist_Binned_PB1(idx,2) = 0;
        idx = find(resTimeHist_Binned_PB1(:,2) <= 0);
        resTimeHist_Binned_PB1_2 = resTimeHist_Binned_PB1;
        resTimeHist_Binned_PB1_2(idx,:) = [];
        % Fit the binned residence time histogram with exponentials;
        [fitpar_PB1,espSigma_PB1, fit1]= ExpDecay_fit_resTime(resTimeHist_Binned_PB1_2, 1 );
        fit_PB1 = [resTimeHist_Binned1(:,1),ExpDecay_fun_resTime(fitpar_PB1,resTimeHist_Binned_PB1(:,1))];

        [fitpar2_PB,espSigma2_PB, fit2]= ExpDecay_2Cmp_fit_resTime(resTimeHist_Binned_PB1_2, [0.1 0.01]);

        Esp_Coef_1 = [fitpar2_PB(1) fitpar2_PB(4)*fitpar2_PB(3)];
        Esp_Fit(:,1) = ExpDecay_fun_resTime(Esp_Coef_1,resTimeHist_Binned_PB1(:,1));

        Esp_Coef_2 = [fitpar2_PB(2) fitpar2_PB(4)*(1 - fitpar2_PB(3))];
        Esp_Fit(:,2) = ExpDecay_fun_resTime(Esp_Coef_2,resTimeHist_Binned_PB1(:,1));

        fit2_PB = [resTimeHist_Binned_PB1(:,1),ExpDecay_2Cmp_fun_resTime(fitpar2_PB,resTimeHist_Binned_PB1(:,1)),Esp_Fit];
        handles.Hist.Res_PB1(:, 1:2) = resTimeHist_Binned_PB1;
        handles.Hist.Res_PB1(:,3) = fit_PB1(:,2);
        handles.Hist.Res_PB1(:,4) = fit2_PB(:,2);
        handles.Hist.Res_PB1(:,5) = fit2_PB(:,3);
        handles.Hist.Res_PB1(:,6) = fit2_PB(:,4);
        [~, pvalRes_PB1,FstatRes_PB1] = FtestModelCompare(resTimeHist_Binned_PB1(:,2),fit_PB1(:,2),fit2_PB(:,2),2,4);
        handles.FitPar.Res_PB1 = [fitpar_PB1, fitpar2_PB, PartialBoundFraction1, pvalRes_PB1; espSigma_PB1, espSigma2_PB, BFerror1, FstatRes_PB1];


        
        axes(handles.C1_C2SurvAxes);
        if strcmp(answer,'Yes');
            resTimeHist_Binned1 = resTimeHist_Binned_PB1;
            fit1 = fit_PB1;
            fit2 = fit2_PB;
            fitpar1 = fitpar_PB1;
            fitpar1_2 = fitpar2_PB;
            espSigma1 = espSigma_PB1;
            espSigma1_2 = espSigma2_PB;
            pvalRes1 = pvalRes_PB1;
            FstatRes1 = FstatRes_PB1;
        end
        if handles.Flag == 1
            bar(resTimeHist_Binned1(:,1), resTimeHist_Binned1(:,2),'EdgeColor', [0 0 0], ...
                'FaceColor',[0 0 0], 'BarWidth', 1, 'LineWidth', 0.75);
            hold on;
            plot(fit1(:,1),fit1(:,2),'g', 'LineWidth', 1);
            plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
            plot(fit2(:,1),fit2(:,3),'b','LineWidth',1);
            plot(fit2(:,1),fit2(:,4),'b--','LineWidth',1);
            hold off;
            xlabel('Residence Time [s]')
            ylabel('Counts');
            title({['Fit of the residence time histogram',...
                ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction1,3), ' \pm ',...
                num2str(BFerror1,2)],...
                ['One Component fit k_{off} = ',num2str(fitpar1(1),3),...
                ' \pm ', num2str(espSigma1(1),2), ' s^{-1}; C_{eq} = ', num2str(fitpar1(2),2)],...
                ['Two Components fit k_{off,1} = ', num2str(fitpar1_2(1),3),...
                ' \pm ', num2str(espSigma1_2(1),2), ' s^{-1}'], ...
                ['k_{off, 2} = ', num2str(fitpar1_2(2),3),...
                ' \pm ', num2str(espSigma1_2(2),2), ' s^{-1}; f_1 = ', num2str(fitpar1_2(3),3),...
                '\pm', num2str(espSigma1_2(3),2), ' C_{eq} = ',  num2str(fitpar1_2(4),2)]});
            legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit','Component 1', 'Component 2');
            box on;
            set(handles.FstatText1,'String',['F stat: ', num2str(FstatRes1,4)]);
            set(handles.pValText1,'String',['p-value: ', num2str(pvalRes1,4)]);

        end




    % end

    % if handles.Flag == 2;
    %     set(handles.What2Plot,'Value',2);
        % Logarithmic sampling of the survival prob
        if handles.Parameters.bin ~= 0

    %         if strcmp(answer,'Yes')
                TimePoints_PB1 = logspace(log10(min(CumTrackLengthHist_PB1(:,1))),...
                    log10(max(CumTrackLengthHist_PB1(:,1))),handles.Parameters.bin);
    %         else
                TimePoints1 = logspace(log10(min(CumTrackLengthHist1(:,1))),...
                    log10(max(CumTrackLengthHist1(:,1))),handles.Parameters.bin);
    %         end
            TrLog1 =[];

            TimePoints1 = round(TimePoints1*1000)/1000;
            TimePoints_PB1 = round(TimePoints_PB1*1000)/1000;
            for i = 1:length(TimePoints1)

                idx = find(CumTrackLengthHist1(:,1) <= TimePoints1(i)+0.0001,1,'last');
                TrLog1(i,:) = CumTrackLengthHist1(idx,:);
            end
            for i = 1:length(TimePoints_PB1)

                idx = find(CumTrackLengthHist_PB1(:,1) <= TimePoints_PB1(i)+0.0001,1,'last');
                TrLog_PB1(i,:) = CumTrackLengthHist_PB1(idx,:);
            end

        else
            TrLog1 = CumTrackLengthHist1;
            TrLog_PB1 = CumTrackLengthHist_PB1;
        end

        tpoint1 = TrLog1(1,1);
        TrLog1(:,1) = TrLog1(:,1) - tpoint1;


    %     ind2fit = 1:length(TrLog);
        % Fit exponentials to the Survival probability
        ind2fit = 1:size(TrLog1,1);
        if handles.Parameters.bin == 0

            [fitpar,espSigma, ~]= ExpDecay_fit(TrLog1(ind2fit,:), 1/(tpoint1*10));
            fit = [TrLog1(:,1),ExpDecay_fun(fitpar,TrLog1(:,1))];
            [fitpar2,espSigma2, ~]= ExpDecay_2Cmp_fit(TrLog1(ind2fit,:), [1/(tpoint1*10) 1/(tpoint1*100)]);
            Esp_Coef_1 = [fitpar2(1) fitpar2(4)*fitpar2(3)];
            Esp_Fit(:,1) = ExpDecay_fun(Esp_Coef_1,TrLog1(:,1));

            Esp_Coef_2 = [fitpar2(2) fitpar2(4)*(1 - fitpar2(3))];
            Esp_Fit(:,2) = ExpDecay_fun(Esp_Coef_2,TrLog1(:,1));
            fit2 = [TrLog1(:,1),ExpDecay_2Cmp_fun(fitpar2,TrLog1(:,1)),Esp_Fit];
        else
            [fitpar,espSigma, fit]= ExpDecay_fit(TrLog1(ind2fit,:), 1/(tpoint1*10));
            [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit(TrLog1(ind2fit,:), [1/(tpoint1*10) 1/(tpoint1*100)]);
        end

        fit(:,1) = fit(:,1) + tpoint1;
        fit2(:,1) = fit2(:,1) + tpoint1;
        TrLog1(:,1) = TrLog1(:,1) + tpoint1;

        tpoint1 = TrLog_PB1(1,1);
        TrLog_PB1(:,1) = TrLog_PB1(:,1) - tpoint1;



        % Fit exponentials to the Survival probability
        if handles.Parameters.bin == 0
            [fitpar_PB,espSigma_PB, ~]= ExpDecay_fit(TrLog_PB1(ind2fit,:), 1/(tpoint1*10));
            fit_PB = [TrLog_PB1(:,1),ExpDecay_fun(fitpar_PB,TrLog_PB1(:,1))];
            [fitpar2_PB,espSigma2_PB, ~]= ExpDecay_2Cmp_fit(TrLog_PB1(ind2fit,:), [1/(tpoint1*10) 1/(tpoint1*100)]);
            Esp_Coef_1 = [fitpar2_PB(1) fitpar2_PB(4)*fitpar2_PB(3)];
            Esp_Fit(:,1) = ExpDecay_fun(Esp_Coef_1,TrLog_PB1(:,1));

            Esp_Coef_2 = [fitpar2_PB(2) fitpar2_PB(4)*(1 - fitpar2_PB(3))];
            Esp_Fit(:,2) = ExpDecay_fun(Esp_Coef_2,TrLog_PB1(:,1));
            fit2_PB = [TrLog_PB1(:,1),ExpDecay_2Cmp_fun(fitpar2_PB,TrLog_PB1(:,1)),Esp_Fit];
        else
            [fitpar_PB,espSigma_PB, fit_PB]= ExpDecay_fit(TrLog_PB1(ind2fit,:), 1/(tpoint1*10));
            [fitpar2_PB,espSigma2_PB, fit2_PB]= ExpDecay_2Cmp_fit(TrLog_PB1(ind2fit,:), [1/(tpoint1*10) 1/(tpoint1*100)]);
        end
        fit_PB(:,1) = fit_PB(:,1) + tpoint1;
        fit2_PB(:,1) = fit2_PB(:,1) + tpoint1;
        TrLog_PB1(:,1) = TrLog_PB1(:,1) + tpoint1;

        handles.Hist.Surv1(:, 1:2) = TrLog1;
        handles.Hist.Surv1(:,3) = fit(:,2);
        handles.Hist.Surv1(:,4) = fit2(:,2);
        handles.Hist.Surv1(:,5) = fit2(:,3);
        handles.Hist.Surv1(:,6) = fit2(:,4);
        [~, pvalSurv, FstatSurv] = FtestModelCompare(TrLog1(:,2),fit(:,2),fit2(:,2),2,4);
        handles.FitPar.Surv1 = [fitpar, fitpar2, PartialBoundFraction1, pvalSurv; espSigma, espSigma2, BFerror1, FstatSurv];

        handles.Hist.Surv_PB1(:, 1:2) = TrLog_PB1;
        handles.Hist.Surv_PB1(:,3) = fit_PB(:,2);
        handles.Hist.Surv_PB1(:,4) = fit2_PB(:,2);
        handles.Hist.Surv_PB1(:,5) = fit2_PB(:,3);
        handles.Hist.Surv_PB1(:,6) = fit2_PB(:,4);
        [~, pvalSurv_PB, FstatSurv_PB] = FtestModelCompare(TrLog_PB1(:,2),fit_PB(:,2),fit2_PB(:,2),2,4);
        handles.FitPar.Surv_PB1 = [fitpar_PB, fitpar2_PB, PartialBoundFraction1, pvalSurv_PB; espSigma_PB, espSigma2_PB, BFerror1, FstatSurv_PB];
        if strcmp(answer,'Yes');
            TrLog1 = TrLog_PB1;
            fit = fit_PB;
            fit2 = fit2_PB;
            fitpar = fitpar_PB;
            fitpar2 = fitpar2_PB;
            espSigma = espSigma_PB;
            espSigma2 = espSigma2_PB;
            pvalSurv = pvalSurv_PB;
            FstatSurv = FstatSurv_PB;
        end


        if handles.Flag == 2
            
            axes(handles.C1_C2SurvAxes);
            plot(TrLog1(:,1), TrLog1(:,2),'ok');
            hold on;
            plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
            plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
            plot(fit2(:,1),fit2(:,3),'b', 'LineWidth', 1);
            plot(fit2(:,1),fit2(:,4),'b--', 'LineWidth', 1);
            hold off;
            xlabel('Time [s]');
            ylabel('Counts');
            title({['Fit of the survival probability',...
                ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction1,3), ' \pm ',...
                num2str(BFerror1,2)],...
                ['One Component fit k_{off} = ',num2str(fitpar(1),3),...
                ' \pm ', num2str(espSigma(1),2), ' s^{-1}; C_{eq} = ', num2str(fitpar(2),2)],...
                ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
                ' \pm ', num2str(espSigma2(1),2), ' s^{-1}'], ...
                ['k_{off, 2} = ', num2str(fitpar2(2),3),...
                ' \pm ', num2str(espSigma2(2),2), ' s^{-1}; f_1 = ', num2str(fitpar2(3),3),...
                '\pm', num2str(espSigma2(3),2), 'C_{eq} = ',  num2str(fitpar2(4),2)]});
            legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit');
            box on;
            set(handles.FstatText1,'String',['F stat: ', num2str(FstatSurv,4)]);
            set(handles.pValText1,'String',['p-value: ', num2str(pvalSurv,4)]);
        end





%     set(handles.totalTracksText1,'String',['Total Tracks: ',num2str(handles.totalTracks1)]);
%     set(handles.boundParticlesText1,'String',['Bound Particles: ',num2str(handles.BoundMolecules1)]);
%     set(handles.totalParticlesText1,'String',['Total Particles: ',num2str(handles.TotalMolecules1)]);
    set(handles.totalROIsText,'String',['Total ROIs: ', num2str(handles.nROIs)]);
    clear Esp_Fit;
    % Remove values < 0 from the histogram for fitting only
        idx1 = find(resTimeHist_Binned2(:,2) < 0);
        resTimeHist_Binned2(idx,2) = 0;
        idx1 = find(resTimeHist_Binned2(:,2) <= 0);
        resTimeHist_Binned1_2 = resTimeHist_Binned2;
        resTimeHist_Binned1_2(idx,:) = [];
        % Fit the binned residence time histogram with exponentials;
        [fitpar1,espSigma1, ~]= ExpDecay_fit_resTime(resTimeHist_Binned1_2, 1 );
        fit1 = [resTimeHist_Binned2(:,1),ExpDecay_fun_resTime(fitpar1,resTimeHist_Binned2(:,1))];
        [fitpar1_2,espSigma1_2, ~]= ExpDecay_2Cmp_fit_resTime(resTimeHist_Binned1_2, [0.1 0.01]);

        Esp_Coef1_1 = [fitpar1_2(1) fitpar1_2(4)*fitpar1_2(3)];
        Esp_Fit1 = [];
        Esp_Fit1(:,1) = ExpDecay_fun_resTime(Esp_Coef1_1,resTimeHist_Binned2(:,1));

        Esp_Coef1_2 = [fitpar1_2(2) fitpar1_2(4)*(1 - fitpar1_2(3))];
        Esp_Fit1(:,2) = ExpDecay_fun_resTime(Esp_Coef1_2,resTimeHist_Binned2(:,1));

        fit2 = [resTimeHist_Binned2(:,1),ExpDecay_2Cmp_fun_resTime(fitpar1_2,resTimeHist_Binned2(:,1)),Esp_Fit1];

        handles.Hist.Res2(:, 1:2) = resTimeHist_Binned2;
        handles.Hist.Res2(:,3) = fit1(:,2);
        handles.Hist.Res2(:,4) = fit2(:,2);
        handles.Hist.Res2(:,5) = fit2(:,3);
        handles.Hist.Res2(:,6) = fit2(:,4);

        [~, pvalRes2,FstatRes2] = FtestModelCompare(resTimeHist_Binned2(:,2),fit1(:,2),fit2(:,2),2,4);

        handles.FitPar.Res2 = [fitpar1, fitpar1_2, PartialBoundFraction2, pvalRes2; espSigma1, espSigma1_2, BFerror2,FstatRes2];
        idx = find(resTimeHist_Binned_PB2(:,2) < 0);
        resTimeHist_Binned_PB2(idx,2) = 0;
        idx = find(resTimeHist_Binned_PB2(:,2) <= 0);
        resTimeHist_Binned_PB1_2 = resTimeHist_Binned_PB2;
        resTimeHist_Binned_PB1_2(idx,:) = [];
        % Fit the binned residence time histogram with exponentials;
        [fitpar_PB1,espSigma_PB1, fit1]= ExpDecay_fit_resTime(resTimeHist_Binned_PB1_2, 1 );
        fit_PB1 = [resTimeHist_Binned_PB2(:,1),ExpDecay_fun_resTime(fitpar_PB1,resTimeHist_Binned_PB2(:,1))];

        [fitpar2_PB,espSigma2_PB, fit2]= ExpDecay_2Cmp_fit_resTime(resTimeHist_Binned_PB1_2, [0.1 0.01]);

        Esp_Coef_1 = [fitpar2_PB(1) fitpar2_PB(4)*fitpar2_PB(3)];
        Esp_Fit(:,1) = ExpDecay_fun_resTime(Esp_Coef_1,resTimeHist_Binned_PB2(:,1));

        Esp_Coef_2 = [fitpar2_PB(2) fitpar2_PB(4)*(1 - fitpar2_PB(3))];
        Esp_Fit(:,2) = ExpDecay_fun_resTime(Esp_Coef_2,resTimeHist_Binned_PB2(:,1));

        fit2_PB = [resTimeHist_Binned_PB2(:,1),ExpDecay_2Cmp_fun_resTime(fitpar2_PB,resTimeHist_Binned_PB2(:,1)),Esp_Fit];
        handles.Hist.Res_PB2(:, 1:2) = resTimeHist_Binned_PB2;
        handles.Hist.Res_PB2(:,3) = fit_PB1(:,2);
        handles.Hist.Res_PB2(:,4) = fit2_PB(:,2);
        handles.Hist.Res_PB2(:,5) = fit2_PB(:,3);
        handles.Hist.Res_PB2(:,6) = fit2_PB(:,4);
        [~, pvalRes_PB2,FstatRes_PB2] = FtestModelCompare(resTimeHist_Binned_PB2(:,2),fit_PB1(:,2),fit2_PB(:,2),2,4);
        handles.FitPar.Res_PB2 = [fitpar_PB1, fitpar2_PB, PartialBoundFraction2, pvalRes_PB2; espSigma_PB1, espSigma2_PB, BFerror2, FstatRes_PB2];


        
        axes(handles.C2_C1SurvAxes);
        if strcmp(answer,'Yes');
            resTimeHist_Binned2 = resTimeHist_Binned_PB2;
            fit1 = fit_PB1;
            fit2 = fit2_PB;
            fitpar1 = fitpar_PB1;
            fitpar1_2 = fitpar2_PB;
            espSigma1 = espSigma_PB1;
            espSigma1_2 = espSigma2_PB;
            pvalRes2 = pvalRes_PB2;
            FstatRes2 = FstatRes_PB2;
        end
        if handles.Flag == 1
            bar(resTimeHist_Binned2(:,1), resTimeHist_Binned2(:,2),'EdgeColor', [0 0 0], ...
                'FaceColor',[0 0 0], 'BarWidth', 1, 'LineWidth', 0.75);
            hold on;
            plot(fit1(:,1),fit1(:,2),'g', 'LineWidth', 1);
            plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
            plot(fit2(:,1),fit2(:,3),'b','LineWidth',1);
            plot(fit2(:,1),fit2(:,4),'b--','LineWidth',1);
            hold off;
            xlabel('Residence Time [s]')
            ylabel('Counts');
            title({['Fit of the residence time histogram',...
                ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction2,3), ' \pm ',...
                num2str(BFerror2,2)],...
                ['One Component fit k_{off} = ',num2str(fitpar1(1),3),...
                ' \pm ', num2str(espSigma1(1),2), ' s^{-1}; C_{eq} = ', num2str(fitpar1(2),2)],...
                ['Two Components fit k_{off,1} = ', num2str(fitpar1_2(1),3),...
                ' \pm ', num2str(espSigma1_2(1),2), ' s^{-1}'], ...
                ['k_{off, 2} = ', num2str(fitpar1_2(2),3),...
                ' \pm ', num2str(espSigma1_2(2),2), ' s^{-1}; f_1 = ', num2str(fitpar1_2(3),3),...
                '\pm', num2str(espSigma1_2(3),2), ' C_{eq} = ',  num2str(fitpar1_2(4),2)]});
            legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit','Component 1', 'Component 2');
            box on;
            set(handles.FstatText2,'String',['F stat: ', num2str(FstatRes2,4)]);
            set(handles.pValText2,'String',['p-value: ', num2str(pvalRes2,4)]);

        end




    % end

    % if handles.Flag == 2;
    %     set(handles.What2Plot,'Value',2);
        % Logarithmic sampling of the survival prob
        clear TrLog1;
        clear TrLog_PB1;
        if handles.Parameters.bin ~= 0

    %         if strcmp(answer,'Yes')
                TimePoints_PB1 = logspace(log10(min(CumTrackLengthHist_PB2(:,1))),...
                    log10(max(CumTrackLengthHist_PB2(:,1))),handles.Parameters.bin);
    %         else
                TimePoints1 = logspace(log10(min(CumTrackLengthHist2(:,1))),...
                    log10(max(CumTrackLengthHist2(:,1))),handles.Parameters.bin);
    %         end
            TrLog1 =[];

            TimePoints1 = round(TimePoints1*1000)/1000;
            TimePoints_PB1 = round(TimePoints_PB1*1000)/1000;
            for i = 1:length(TimePoints1)

                idx = find(CumTrackLengthHist2(:,1) <= TimePoints1(i)+0.0001,1,'last');
                TrLog1(i,:) = CumTrackLengthHist2(idx,:);
            end
            for i = 1:length(TimePoints_PB1)

                idx = find(CumTrackLengthHist_PB2(:,1) <= TimePoints_PB1(i)+0.0001,1,'last');
                TrLog_PB1(i,:) = CumTrackLengthHist_PB2(idx,:);
            end

        else
            TrLog1 = CumTrackLengthHist2;
            TrLog_PB1 = CumTrackLengthHist_PB2;
        end

        tpoint1 = TrLog1(1,1);
        TrLog1(:,1) = TrLog1(:,1) - tpoint1;


    %     ind2fit = 1:length(TrLog);
        % Fit exponentials to the Survival probability
        ind2fit = 1:size(TrLog1,1);
        if handles.Parameters.bin == 0

            [fitpar,espSigma, ~]= ExpDecay_fit(TrLog1(ind2fit,:), 1/(tpoint1*10));
            fit = [TrLog1(:,1),ExpDecay_fun(fitpar,TrLog1(:,1))];
            [fitpar2,espSigma2, ~]= ExpDecay_2Cmp_fit(TrLog1(ind2fit,:), [1/(tpoint1*10) 1/(tpoint1*100)]);
            Esp_Coef_1 = [fitpar2(1) fitpar2(4)*fitpar2(3)];
            Esp_Fit(:,1) = ExpDecay_fun(Esp_Coef_1,TrLog1(:,1));

            Esp_Coef_2 = [fitpar2(2) fitpar2(4)*(1 - fitpar2(3))];
            Esp_Fit(:,2) = ExpDecay_fun(Esp_Coef_2,TrLog1(:,1));
            fit2 = [TrLog1(:,1),ExpDecay_2Cmp_fun(fitpar2,TrLog1(:,1)),Esp_Fit];
        else
            [fitpar,espSigma, fit]= ExpDecay_fit(TrLog1(ind2fit,:), 1/(tpoint1*10));
            [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit(TrLog1(ind2fit,:), [1/(tpoint1*10) 1/(tpoint1*100)]);
        end

        fit(:,1) = fit(:,1) + tpoint1;
        fit2(:,1) = fit2(:,1) + tpoint1;
        TrLog1(:,1) = TrLog1(:,1) + tpoint1;

        tpoint1 = TrLog_PB1(1,1);
        TrLog_PB1(:,1) = TrLog_PB1(:,1) - tpoint1;



        % Fit exponentials to the Survival probability
        if handles.Parameters.bin == 0
            [fitpar_PB,espSigma_PB, ~]= ExpDecay_fit(TrLog_PB1(ind2fit,:), 1/(tpoint1*10));
            fit_PB = [TrLog_PB1(:,1),ExpDecay_fun(fitpar_PB,TrLog_PB1(:,1))];
            [fitpar2_PB,espSigma2_PB, ~]= ExpDecay_2Cmp_fit(TrLog_PB1(ind2fit,:), [1/(tpoint1*10) 1/(tpoint1*100)]);
            Esp_Coef_1 = [fitpar2_PB(1) fitpar2_PB(4)*fitpar2_PB(3)];
            Esp_Fit(:,1) = ExpDecay_fun(Esp_Coef_1,TrLog_PB1(:,1));

            Esp_Coef_2 = [fitpar2_PB(2) fitpar2_PB(4)*(1 - fitpar2_PB(3))];
            Esp_Fit(:,2) = ExpDecay_fun(Esp_Coef_2,TrLog_PB1(:,1));
            fit2_PB = [TrLog_PB1(:,1),ExpDecay_2Cmp_fun(fitpar2_PB,TrLog_PB1(:,1)),Esp_Fit];
        else
            [fitpar_PB,espSigma_PB, fit_PB]= ExpDecay_fit(TrLog_PB1(ind2fit,:), 1/(tpoint1*10));
            [fitpar2_PB,espSigma2_PB, fit2_PB]= ExpDecay_2Cmp_fit(TrLog_PB1(ind2fit,:), [1/(tpoint1*10) 1/(tpoint1*100)]);
        end
        fit_PB(:,1) = fit_PB(:,1) + tpoint1;
        fit2_PB(:,1) = fit2_PB(:,1) + tpoint1;
        TrLog_PB1(:,1) = TrLog_PB1(:,1) + tpoint1;

        handles.Hist.Surv2(:, 1:2) = TrLog1;
        handles.Hist.Surv2(:,3) = fit(:,2);
        handles.Hist.Surv2(:,4) = fit2(:,2);
        handles.Hist.Surv2(:,5) = fit2(:,3);
        handles.Hist.Surv2(:,6) = fit2(:,4);
        [~, pvalSurv2, FstatSurv2] = FtestModelCompare(TrLog1(:,2),fit(:,2),fit2(:,2),2,4);
        handles.FitPar.Surv2 = [fitpar, fitpar2, PartialBoundFraction1, pvalSurv2; espSigma, espSigma2, BFerror1, FstatSurv2];

        handles.Hist.Surv_PB2(:, 1:2) = TrLog_PB1;
        handles.Hist.Surv_PB2(:,3) = fit_PB(:,2);
        handles.Hist.Surv_PB2(:,4) = fit2_PB(:,2);
        handles.Hist.Surv_PB2(:,5) = fit2_PB(:,3);
        handles.Hist.Surv_PB2(:,6) = fit2_PB(:,4);
        [~, pvalSurv_PB2, FstatSurv_PB2] = FtestModelCompare(TrLog_PB1(:,2),fit_PB(:,2),fit2_PB(:,2),2,4);
        handles.FitPar.Surv_PB2 = [fitpar_PB, fitpar2_PB, PartialBoundFraction2, pvalSurv_PB2; espSigma_PB, espSigma2_PB, BFerror2, FstatSurv_PB2];
        if strcmp(answer,'Yes');
            TrLog1 = TrLog_PB1;
            fit = fit_PB;
            fit2 = fit2_PB;
            fitpar = fitpar_PB;
            fitpar2 = fitpar2_PB;
            espSigma = espSigma_PB;
            espSigma2 = espSigma2_PB;
            pvalSurv2 = pvalSurv_PB2;
            FstatSurv2 = FstatSurv_PB2;
        end


        if handles.Flag == 2
            
            axes(handles.C2_C1SurvAxes);
            plot(TrLog1(:,1), TrLog1(:,2),'ok');
            hold on;
            plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
            plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
            plot(fit2(:,1),fit2(:,3),'b', 'LineWidth', 1);
            plot(fit2(:,1),fit2(:,4),'b--', 'LineWidth', 1);
            hold off;
            xlabel('Time [s]');
            ylabel('Counts');
            title({['Fit of the survival probability',...
                ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction2,3), ' \pm ',...
                num2str(BFerror2,2)],...
                ['One Component fit k_{off} = ',num2str(fitpar(1),3),...
                ' \pm ', num2str(espSigma(1),2), ' s^{-1}; C_{eq} = ', num2str(fitpar(2),2)],...
                ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
                ' \pm ', num2str(espSigma2(1),2), ' s^{-1}'], ...
                ['k_{off, 2} = ', num2str(fitpar2(2),3),...
                ' \pm ', num2str(espSigma2(2),2), ' s^{-1}; f_1 = ', num2str(fitpar2(3),3),...
                '\pm', num2str(espSigma2(3),2), 'C_{eq} = ',  num2str(fitpar2(4),2)]});
            legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit');
            box on;
            set(handles.FstatText2,'String',['F stat: ', num2str(FstatSurv2,4)]);
            set(handles.pValText2,'String',['p-value: ', num2str(pvalSurv2,4)]);
        end





    set(handles.totalTracksText,'String',['Total Tracks: Ch1: ',num2str(handles.totalTracks1), ', Ch2: ',num2str(handles.totalTracks2) ]);
    set(handles.boundParticlesText,'String',['Bound Particles: Ch1: ',num2str(handles.BoundMolecules1), ', Ch2: ',num2str(handles.BoundMolecules2)]);
    set(handles.totalParticlesText,'String',['Total Particles: Ch1: ',num2str(handles.TotalMolecules1), ' Ch2: ',num2str(handles.TotalMolecules2)]);
    set(handles.totalColocPartText,'String',['Total Coloc Particles: ', num2str(handles.ColocParticles)]);
    set(handles.totalROIsText,'String',['Total ROIs: ', num2str(handles.nROIs)]);

    figname = get(gcf,'Name');
    figname = [figname ' ' ROIChoice{1}];
    set(gcf,'Name',figname);
%     set(handles.PieChart,'Enable','on');
  
    set(handles.frame2frameMaxText,'String',['Max. Frame to Frame Jump: ', num2str(handles.Parameters.ThreshL), ' um']);
    set(handles.End2EndMaxText,'String',['Max. End to End Jump: ', num2str(handles.Parameters.ThreshH), ' um']);
    set(handles.NminText,'String',['Min. Bound Frames (Nmin): ', num2str(handles.Parameters.minBoundFrames)]);
    set(handles.MinInteractText,'String',['Min. Interaction Frames: ', num2str(handles.Parameters.minInteract)]);
    set(handles.MaxSep2colorText,'String',['Max. Color Separation: ', num2str(handles.Parameters.Dist2color)]);
    
    
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TwoColorBinding_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TwoColorBinding_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = [];


% --- Executes on selection change in What2Plot.
function What2Plot_Callback(hObject, eventdata, handles)
% hObject    handle to What2Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns What2Plot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from What2Plot
PartialBoundFraction1 = handles.PartBoundFrac1;
BFerror1 = handles.BFerror1;
if get(hObject,'Value') == 1;
    

    if get(handles.PBCorrect,'Value') == 1 
        Hist1 = handles.Hist.Res_PB1;
        pval1 = handles.FitPar.Res_PB1(1,8);
        Fstat1 = handles.FitPar.Res_PB1(2,8);
        
    else
        Hist1 = handles.Hist.Res1;
        pval1 = handles.FitPar.Res1(1,8);
        Fstat1 = handles.FitPar.Res1(2,8);
    end
    
    if size(Hist1,2) < 3

        % Bin the histogram according to the selected settings;

        binFactor = handles.Parameters.bin + 1;
        nBins1 = floor(length(Hist1(:,1))/binFactor);
        
            TimeBins1 = reshape(Hist1(1:nBins1*binFactor,1), binFactor, nBins1);
            CountsBins1 = reshape(Hist1(1:nBins1*binFactor,2), binFactor, nBins1);
        

        if binFactor ~= 1
            resTimeHist_Binned1 = [];
            resTimeHist_Binned1(:,1) = mean(TimeBins1);
            resTimeHist_Binned1(:,2) = sum(CountsBins1);
        else
            resTimeHist_Binned1 = Hist1;
        end

        % Remove values < 0 from the histogram
        idx = find(resTimeHist_Binned1(:,2) <= 0);
        resTimeHist_Binned1(idx,:) = [];
        % Fit the binned residence time histogram with exponentials;
        [fitpar,espSigma, fit]= ExpDecay_fit_resTime(resTimeHist_Binned1, 1 );
        [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit_resTime(resTimeHist_Binned1, [1 0.1]);
        if get(handles.PBCorrect,'Value') == 1 
            handles.Hist.Res_PB1(:, 1:2) = resTimeHist_Binned1;
            handles.Hist.Res_PB1(:,3) = fit(:,2);
            handles.Hist.Res_PB1(:,4) = fit2(:,2);
            handles.Hist.Res_PB1(:,5) = fit2(:,3);
            handles.Hist.Res_PB1(:,6) = fit2(:,4);
            
            handles.FitPar.Res_PB1 = [fitpar, fitpar2, PartialBoundFraction1; espSigma, espSigma2, BFerror1];
        else
            handles.Hist.Res1(:, 1:2) = resTimeHist_Binned1;
            handles.Hist.Res1(:,3) = fit(:,2);
            handles.Hist.Res1(:,4) = fit2(:,2);
            handles.Hist.Res1(:,5) = fit2(:,3);
            handles.Hist.Res1(:,6) = fit2(:,4);
            
            handles.FitPar.Res1 = [fitpar, fitpar2, PartialBoundFraction1; espSigma, espSigma2, BFerror1];
        end
        
    else
        resTimeHist_Binned1 = Hist1;
        fit = [Hist1(:,1), Hist1(:,3)];
        fit2 = [Hist1(:,1),Hist1(:,4:6)];
        if get(handles.PBCorrect,'Value') == 1 
            fitpar = handles.FitPar.Res_PB1(1,1:2);
            espSigma = handles.FitPar.Res_PB1(2,1:2);
            fitpar2 = handles.FitPar.Res_PB1(1,3:6);
            espSigma2 = handles.FitPar.Res_PB1(2,3:6);
        else
            fitpar = handles.FitPar.Res1(1,1:2);
            espSigma = handles.FitPar.Res1(2,1:2);
            fitpar2 = handles.FitPar.Res1(1,3:6);
            espSigma2 = handles.FitPar.Res1(2,3:6);
        end
    end
        
    
    
    
%     HistAxes = findobj ('Tag', 'hist_axes');
    axes(handles.C1_C2SurvAxes);
    
    bar(resTimeHist_Binned1(:,1), resTimeHist_Binned1(:,2),'EdgeColor', [0 0 0], ...
        'FaceColor',[0 0 0], 'BarWidth', 1, 'LineWidth', 0.75);
    hold on;
    plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,3),'b','LineWidth',1);
    plot(fit2(:,1),fit2(:,4),'b--','LineWidth',1);
    hold off;
    xlabel('Residence Time [s]')
    ylabel('Counts');
    title({['Fit of the residence time histogram',...
        ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction1,3), ' \pm ',...
        num2str(BFerror1,2)],...
        ['One Component fit k_{off} = ',num2str(fitpar(1),3),...
        ' \pm ', num2str(espSigma(1),2), ' s^{-1}; C_{eq} = ', num2str(fitpar(2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
        ' \pm ', num2str(espSigma2(1),2), ' s^{-1}'], ...
        ['k_{off, 2} = ', num2str(fitpar2(2),3),...
        ' \pm ', num2str(espSigma2(2),2), ' s^{-1}; f_1 = ', num2str(fitpar2(3),3),...
        '\pm', num2str(espSigma2(3),2), ' C_{eq} = ',  num2str(fitpar2(4),2)]});
    legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit','Component 1', 'Component 2');
    box on;
    
    
    
    
    


else
    if get(handles.PBCorrect,'Value') == 1 
        Hist1 = handles.Hist.Surv_PB1;
        pval1 = handles.FitPar.Surv_PB1(1,8);
        Fstat1 = handles.FitPar.Surv_PB1(2,8);
    else
        Hist1 = handles.Hist.Surv1;
        pval1 = handles.FitPar.Surv1(1,8);
        Fstat1 = handles.FitPar.Surv1(2,8);
    end
    
    if size(Hist1,2) < 3
    
        % Logarithmic sampling of the survival prob
        if handles.Parameters.bin ~= 0
            
            
            TimePoints1 = logspace(log10(min(Hist1(:,1))),...
                log10(max(Hist1(:,1))),handles.Parameters.bin);
            
            TrLog1 =[];
            
            TimePoints1 = round(TimePoints1*1000)/1000;
            for i = 1:length(TimePoints1)
                
                idx = find(Hist1(:,1) <= TimePoints1(i)+0.0001,1,'last');
                TrLog1(i,:) = Hist1(idx,:);
            end
            
        else
            TrLog1 = Hist1;
        end
        
        tpoint1 = TrLog1(1,1);
        TrLog1(:,1) = TrLog1(:,1) - tpoint1;
        
        
        
        
        % Fit exponentials to the Survival probability
        [fitpar,espSigma, fit]= ExpDecay_fit(TrLog1, 1/(tpoint1*10));
        [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit(TrLog1, [1/(tpoint1*10) 1/(tpoint1*100)]);
        fit(:,1) = fit(:,1) + tpoint1;
        fit2(:,1) = fit2(:,1) + tpoint1;
        TrLog1(:,1) = TrLog1(:,1) + tpoint1;
        if get(handles.PBCorrect,'Value') == 0
            handles.Hist.Surv1(:, 1:2) = TrLog1;
            handles.Hist.Surv1(:,3) = fit(:,2);
            handles.Hist.Surv1(:,4) = fit2(:,2);
            handles.Hist.Surv1(:,5) = fit2(:,3);
            handles.Hist.Surv1(:,6) = fit2(:,4);
            handles.FitPar.Surv1 = [fitpar, fitpar2, PartialBoundFraction1; espSigma, espSigma2, BFerror1];
        else
            handles.Hist.Surv_PB1(:, 1:2) = TrLog1;
            handles.Hist.Surv_PB1(:,3) = fit(:,2);
            handles.Hist.Surv_PB1(:,4) = fit2(:,2);
            handles.Hist.Surv_PB1(:,5) = fit2(:,3);
            handles.Hist.Surv_PB1(:,6) = fit2(:,4);
            handles.FitPar.Surv_PB1 = [fitpar, fitpar2, PartialBoundFraction1; espSigma, espSigma2, BFerror1];
        end
    else
        TrLog1 = Hist1;
        fit = [Hist1(:,1),Hist1(:,3)];
        fit2 = [Hist1(:,1),Hist1(:,4:6)];
        if get(handles.PBCorrect,'Value') == 1
            fitpar = handles.FitPar.Surv_PB1(1,1:2);
            espSigma = handles.FitPar.Surv_PB1(2,1:2);
            fitpar2 = handles.FitPar.Surv_PB1(1,3:6);
            espSigma2 = handles.FitPar.Surv_PB1(2,3:6);
        else
            fitpar = handles.FitPar.Surv1(1,1:2);
            espSigma = handles.FitPar.Surv1(2,1:2);
            fitpar2 = handles.FitPar.Surv1(1,3:6);
            espSigma2 = handles.FitPar.Surv1(2,3:6);
        end
    end
    
%     HistAxes = findobj ('Tag', '');
    axes(handles.C1_C2SurvAxes);
    plot(TrLog1(:,1), TrLog1(:,2),'ok');
    hold on;
    plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,3),'b', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,4),'b--', 'LineWidth', 1);
    hold off;
    xlabel('Time [s]');
    ylabel('Counts');
    title({['Fit of the survival probability',...
        ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction1,3), ' \pm ',...
        num2str(BFerror1,2)],...
        ['One Component fit k_{off} = ',num2str(fitpar(1),3),...
        ' \pm ', num2str(espSigma(1),2), ' s^{-1}; C_{eq} = ', num2str(fitpar(2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
        ' \pm ', num2str(espSigma2(1),2), ' s^{-1}'], ...
        ['k_{off, 2} = ', num2str(fitpar2(2),3),...
        ' \pm ', num2str(espSigma2(2),2), ' s^{-1}; f_1 = ', num2str(fitpar2(3),3),...
        '\pm', num2str(espSigma2(3),2), 'C_{eq} = ',  num2str(fitpar2(4),2)]});
    legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit');
    box on;
    
    
        
    
    
end
handles.CumHist1 = Hist1;
set(handles.FstatText1,'String',['F stat: ', num2str(Fstat1,4)]);
set(handles.pValText1,'String',['p-value: ', num2str(pval1,4)]);


%Channel 2
PartialBoundFraction2 = handles.PartBoundFrac2;
BFerror2 = handles.BFerror2;
if get(hObject,'Value') == 1;
    

    if get(handles.PBCorrect,'Value') == 1 
        Hist2 = handles.Hist.Res_PB2;
        pval2 = handles.FitPar.Res_PB2(1,8);
        Fstat2 = handles.FitPar.Res_PB2(2,8);
        
    else
        Hist2 = handles.Hist.Res2;
        pval2 = handles.FitPar.Res2(1,8);
        Fstat2 = handles.FitPar.Res2(2,8);
    end
    
    if size(Hist2,2) < 3

        % Bin the histogram according to the selected settings;

        binFactor = handles.Parameters.bin + 1;
        nBins2 = floor(length(Hist2(:,1))/binFactor);
        
            TimeBins2 = reshape(Hist2(1:nBins2*binFactor,1), binFactor, nBins2);
            CountsBins2 = reshape(Hist2(1:nBins2*binFactor,2), binFactor, nBins2);
        

        if binFactor ~= 1
            resTimeHist_Binned2 = [];
            resTimeHist_Binned2(:,1) = mean(TimeBins2);
            resTimeHist_Binned2(:,2) = sum(CountsBins2);
        else
            resTimeHist_Binned2 = Hist2;
        end

        % Remove values < 0 from the histogram
        idx = find(resTimeHist_Binned2(:,2) <= 0);
        resTimeHist_Binned2(idx,:) = [];
        % Fit the binned residence time histogram with exponentials;
        [fitpar,espSigma, fit]= ExpDecay_fit_resTime(resTimeHist_Binned2, 1 );
        [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit_resTime(resTimeHist_Binned2, [1 0.1]);
        if get(handles.PBCorrect,'Value') == 1 
            handles.Hist.Res_PB2(:, 1:2) = resTimeHist_Binned2;
            handles.Hist.Res_PB2(:,3) = fit(:,2);
            handles.Hist.Res_PB2(:,4) = fit2(:,2);
            handles.Hist.Res_PB2(:,5) = fit2(:,3);
            handles.Hist.Res_PB2(:,6) = fit2(:,4);
            
            handles.FitPar.Res_PB2 = [fitpar, fitpar2, PartialBoundFraction2; espSigma, espSigma2, BFerror2];
        else
            handles.Hist.Res2(:, 1:2) = resTimeHist_Binned2;
            handles.Hist.Res2(:,3) = fit(:,2);
            handles.Hist.Res2(:,4) = fit2(:,2);
            handles.Hist.Res2(:,5) = fit2(:,3);
            handles.Hist.Res2(:,6) = fit2(:,4);
            
            handles.FitPar.Res2 = [fitpar, fitpar2, PartialBoundFraction2; espSigma, espSigma2, BFerror2];
        end
        
    else
        resTimeHist_Binned2 = Hist2;
        fit = [Hist2(:,1), Hist2(:,3)];
        fit2 = [Hist2(:,1),Hist2(:,4:6)];
        if get(handles.PBCorrect,'Value') == 1 
            fitpar = handles.FitPar.Res_PB2(1,1:2);
            espSigma = handles.FitPar.Res_PB2(2,1:2);
            fitpar2 = handles.FitPar.Res_PB2(1,3:6);
            espSigma2 = handles.FitPar.Res_PB2(2,3:6);
        else
            fitpar = handles.FitPar.Res2(1,1:2);
            espSigma = handles.FitPar.Res2(2,1:2);
            fitpar2 = handles.FitPar.Res2(1,3:6);
            espSigma2 = handles.FitPar.Res2(2,3:6);
        end
    end
        
    
    
    
%     HistAxes = findobj ('Tag', 'hist_axes');
    axes(handles.C2_C1SurvAxes);
    
    bar(resTimeHist_Binned2(:,1), resTimeHist_Binned2(:,2),'EdgeColor', [0 0 0], ...
        'FaceColor',[0 0 0], 'BarWidth', 1, 'LineWidth', 0.75);
    hold on;
    plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,3),'b','LineWidth',1);
    plot(fit2(:,1),fit2(:,4),'b--','LineWidth',1);
    hold off;
    xlabel('Residence Time [s]')
    ylabel('Counts');
    title({['Fit of the residence time histogram',...
        ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction2,3), ' \pm ',...
        num2str(BFerror2,2)],...
        ['One Component fit k_{off} = ',num2str(fitpar(1),3),...
        ' \pm ', num2str(espSigma(1),2), ' s^{-1}; C_{eq} = ', num2str(fitpar(2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
        ' \pm ', num2str(espSigma2(1),2), ' s^{-1}'], ...
        ['k_{off, 2} = ', num2str(fitpar2(2),3),...
        ' \pm ', num2str(espSigma2(2),2), ' s^{-1}; f_1 = ', num2str(fitpar2(3),3),...
        '\pm', num2str(espSigma2(3),2), ' C_{eq} = ',  num2str(fitpar2(4),2)]});
    legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit','Component 1', 'Component 2');
    box on;
    
    
    
    
    


else
    if get(handles.PBCorrect,'Value') == 1 
        Hist2 = handles.Hist.Surv_PB2;
        pval2 = handles.FitPar.Surv_PB2(1,8);
        Fstat2 = handles.FitPar.Surv_PB2(2,8);
    else
        Hist2 = handles.Hist.Surv2;
        pval2 = handles.FitPar.Surv2(1,8);
        Fstat2 = handles.FitPar.Surv2(2,8);
    end
    
    if size(Hist2,2) < 3
    
        % Logarithmic sampling of the survival prob
        if handles.Parameters.bin ~= 0
            
            
            TimePoints2 = logspace(log10(min(Hist2(:,1))),...
                log10(max(Hist2(:,1))),handles.Parameters.bin);
            
            TrLog2 =[];
            
            TimePoints2 = round(TimePoints2*1000)/1000;
            for i = 1:length(TimePoints2)
                
                idx = find(Hist2(:,1) <= TimePoints2(i)+0.0001,1,'last');
                TrLog2(i,:) = Hist2(idx,:);
            end
            
        else
            TrLog2 = Hist2;
        end
        
        tpoint1 = TrLog2(1,1);
        TrLog2(:,1) = TrLog2(:,1) - tpoint1;
        
        
        
        
        % Fit exponentials to the Survival probability
        [fitpar,espSigma, fit]= ExpDecay_fit(TrLog2, 1/(tpoint1*10));
        [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit(TrLog2, [1/(tpoint1*10) 1/(tpoint1*100)]);
        fit(:,1) = fit(:,1) + tpoint1;
        fit2(:,1) = fit2(:,1) + tpoint1;
        TrLog2(:,1) = TrLog2(:,1) + tpoint1;
        if get(handles.PBCorrect,'Value') == 0
            handles.Hist.Surv2(:, 1:2) = TrLog1;
            handles.Hist.Surv2(:,3) = fit(:,2);
            handles.Hist.Surv2(:,4) = fit2(:,2);
            handles.Hist.Surv2(:,5) = fit2(:,3);
            handles.Hist.Surv2(:,6) = fit2(:,4);
            handles.FitPar.Surv2 = [fitpar, fitpar2, PartialBoundFraction2; espSigma, espSigma2, BFerror2];
        else
            handles.Hist.Surv_PB2(:, 1:2) = TrLog2;
            handles.Hist.Surv_PB2(:,3) = fit(:,2);
            handles.Hist.Surv_PB2(:,4) = fit2(:,2);
            handles.Hist.Surv_PB2(:,5) = fit2(:,3);
            handles.Hist.Surv_PB2(:,6) = fit2(:,4);
            handles.FitPar.Surv_PB2 = [fitpar, fitpar2, PartialBoundFraction2; espSigma, espSigma2, BFerror2];
        end
    else
        TrLog2 = Hist2;
        fit = [Hist2(:,1),Hist2(:,3)];
        fit2 = [Hist2(:,1),Hist2(:,4:6)];
        if get(handles.PBCorrect,'Value') == 1
            fitpar = handles.FitPar.Surv_PB2(1,1:2);
            espSigma = handles.FitPar.Surv_PB2(2,1:2);
            fitpar2 = handles.FitPar.Surv_PB2(1,3:6);
            espSigma2 = handles.FitPar.Surv_PB2(2,3:6);
        else
            fitpar = handles.FitPar.Surv2(1,1:2);
            espSigma = handles.FitPar.Surv2(2,1:2);
            fitpar2 = handles.FitPar.Surv2(1,3:6);
            espSigma2 = handles.FitPar.Surv2(2,3:6);
        end
    end
    
%     HistAxes = findobj ('Tag', '');
    axes(handles.C2_C1SurvAxes);
    plot(TrLog2(:,1), TrLog2(:,2),'ok');
    hold on;
    plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,3),'b', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,4),'b--', 'LineWidth', 1);
    hold off;
    xlabel('Time [s]');
    ylabel('Counts');
    title({['Fit of the survival probability',...
        ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction2,3), ' \pm ',...
        num2str(BFerror2,2)],...
        ['One Component fit k_{off} = ',num2str(fitpar(1),3),...
        ' \pm ', num2str(espSigma(1),2), ' s^{-1}; C_{eq} = ', num2str(fitpar(2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
        ' \pm ', num2str(espSigma2(1),2), ' s^{-1}'], ...
        ['k_{off, 2} = ', num2str(fitpar2(2),3),...
        ' \pm ', num2str(espSigma2(2),2), ' s^{-1}; f_1 = ', num2str(fitpar2(3),3),...
        '\pm', num2str(espSigma2(3),2), 'C_{eq} = ',  num2str(fitpar2(4),2)]});
    legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit');
    box on;
    
    
        
    
    
end
handles.CumHist2 = Hist2;
set(handles.FstatText2,'String',['F stat: ', num2str(Fstat2,4)]);
set(handles.pValText2,'String',['p-value: ', num2str(pval2,4)]);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function What2Plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to What2Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PBCorrect.
function PBCorrect_Callback(hObject, eventdata, handles)
% hObject    handle to PBCorrect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PBCorrect
PartialBoundFraction1 = handles.PartBoundFrac1;
BFerror1 = handles.BFerror1;
if get(handles.What2Plot,'Value') == 1;
    

    if get(handles.PBCorrect,'Value') == 1 
        Hist1 = handles.Hist.Res_PB1;
        pval1 = handles.FitPar.Res1(1,8);
        Fstat1 = handles.FitPar.Res1(2,8);
    else
        Hist1 = handles.Hist.Res1;
        pval1 = handles.FitPar.Res_PB1(1,8);
        Fstat1 = handles.FitPar.Res_PB1(2,8);
    end
    
    if size(Hist1,2) < 3

        % Bin the histogram according to the selected settings;

        binFactor = handles.Parameters.bin + 1;
        nBins1 = floor(length(Hist1(:,1))/binFactor);
        
            TimeBins1 = reshape(Hist1(1:nBins1*binFactor,1), binFactor, nBins1);
            CountsBins1 = reshape(Hist1(1:nBins1*binFactor,2), binFactor, nBins1);
        

        if binFactor ~= 1
            resTimeHist_Binned1 = [];
            resTimeHist_Binned1(:,1) = mean(TimeBins1);
            resTimeHist_Binned1(:,2) = sum(CountsBins1);
        else
            resTimeHist_Binned1 = Hist1;
        end

        % Remove values < 0 from the histogram
        idx = find(resTimeHist_Binned1(:,2) <= 0);
        resTimeHist_Binned1(idx,:) = [];
        % Fit the binned residence time histogram with exponentials;
        [fitpar,espSigma, fit]= ExpDecay_fit_resTime(resTimeHist_Binned1, 1 );
        [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit_resTime(resTimeHist_Binned1, [1 0.1]);
        if get(handles.PBCorrect,'Value') == 1 
            handles.Hist.Res_PB1(:, 1:2) = resTimeHist_Binned1;
            handles.Hist.Res_PB1(:,3) = fit(:,2);
            handles.Hist.Res_PB1(:,4) = fit2(:,2);
            handles.Hist.Res_PB1(:,5) = fit2(:,3);
            handles.Hist.Res_PB1(:,6) = fit2(:,4);
            
            handles.FitPar.Res_PB1 = [fitpar, fitpar2, PartialBoundFraction1; espSigma, espSigma2, BFerror1];
        else
            handles.Hist.Res1(:, 1:2) = resTimeHist_Binned1;
            handles.Hist.Res1(:,3) = fit(:,2);
            handles.Hist.Res1(:,4) = fit2(:,2);
            handles.Hist.Res1(:,5) = fit2(:,3);
            handles.Hist.Res1(:,6) = fit2(:,4);
            
            handles.FitPar.Res1 = [fitpar, fitpar2, PartialBoundFraction1; espSigma, espSigma2, BFerror1];
        end
        
    else
        resTimeHist_Binned1 = Hist1;
        fit = [Hist1(:,1),Hist1(:,3)];
        fit2 = [Hist1(:,1),Hist1(:,4:6)];
        if get(handles.PBCorrect,'Value') == 1 
            fitpar = handles.FitPar.Res_PB1(1,1:2);
            espSigma = handles.FitPar.Res_PB1(2,1:2);
            fitpar2 = handles.FitPar.Res_PB1(1,3:6);
            espSigma2 = handles.FitPar.Res_PB1(2,3:6);
        else
            fitpar = handles.FitPar.Res1(1,1:2);
            espSigma = handles.FitPar.Res1(2,1:2);
            fitpar2 = handles.FitPar.Res1(1,3:6);
            espSigma2 = handles.FitPar.Res1(2,3:6);
        end
    end
        
    
    
    
%     HistAxes = findobj ('Tag', 'hist_axes');
    axes(handles.C1_C2SurvAxes);
    bar(resTimeHist_Binned1(:,1), resTimeHist_Binned1(:,2),'EdgeColor', [0 0 0], ...
        'FaceColor',[0 0 0], 'BarWidth', 1, 'LineWidth', 0.75);
    hold on;
    plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,3),'b','LineWidth',1);
    plot(fit2(:,1),fit2(:,4),'b--','LineWidth',1);
    hold off
    xlabel('Residence Time [s]')
    ylabel('Counts');
    title({['Fit of the residence time histogram',...
        ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction1,3), ' \pm ',...
        num2str(BFerror1,2)],...
        ['One Component fit k_{off} = ',num2str(fitpar(1),3),...
        ' \pm ', num2str(espSigma(1),2), ' s^{-1}; C_{eq} = ', num2str(fitpar(2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
        ' \pm ', num2str(espSigma2(1),2), ' s^{-1}'], ...
        ['k_{off, 2} = ', num2str(fitpar2(2),3),...
        ' \pm ', num2str(espSigma2(2),2), ' s^{-1}; f_1 = ', num2str(fitpar2(3),3),...
        '\pm', num2str(espSigma2(3),2), ' C_{eq} = ',  num2str(fitpar2(4),2)]});
    legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit','Component 1', 'Component 2');
    box on;
    
    
    
    
    


else
    if get(handles.PBCorrect,'Value') == 1 
        Hist1 = handles.Hist.Surv_PB1;
        pval1 = handles.FitPar.Surv_PB1(1,8);
        Fstat1 = handles.FitPar.Surv_PB1(2,8);
    else
        Hist1 = handles.Hist.Surv1;
        pval1 = handles.FitPar.Surv1(1,8);
        Fstat1 = handles.FitPar.Surv1(2,8);
    end
    
    if size(Hist1,2) < 3
    
        % Logarithmic sampling of the survival prob
        if handles.Parameters.bin ~= 0
            
            
            TimePoints1 = logspace(log10(min(Hist1(:,1))),...
                log10(max(Hist1(:,1))),handles.Parameters.bin);
            
            TrLog1 =[];
            
            TimePoints1 = round(TimePoints1*1000)/1000;
            for i = 1:length(TimePoints1)
                
                idx = find(Hist1(:,1) <= TimePoints1(i)+0.0001,1,'last');
                TrLog1(i,:) = Hist1(idx,:);
            end
            
        else
            TrLog1 = Hist1;
        end
        
        tpoint1 = TrLog1(1,1);
        TrLog1(:,1) = TrLog1(:,1) - tpoint1;
        
        
        
        
        % Fit exponentials to the Survival probability
        [fitpar,espSigma, fit]= ExpDecay_fit(TrLog1, 1/(tpoint1*10));
        [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit(TrLog1, [1/(tpoint1*10) 1/(tpoint1*100)]);
        fit(:,1) = fit(:,1) + tpoint1;
        fit2(:,1) = fit2(:,1) + tpoint1;
        TrLog1(:,1) = TrLog1(:,1) + tpoint1;
        if get(handles.PBCorrect,'Value') == 0
            handles.Hist.Surv1(:, 1:2) = TrLog1;
            handles.Hist.Surv1(:,3) = fit(:,2);
            handles.Hist.Surv1(:,4) = fit2(:,2);
            handles.Hist.Surv1(:,5) = fit2(:,3);
            handles.Hist.Surv1(:,6) = fit2(:,4);
            handles.FitPar.Surv1 = [fitpar, fitpar2, PartialBoundFraction1; espSigma, espSigma2, BFerror1];
        else
            handles.Hist.Surv_PB1(:, 1:2) = TrLog1;
            handles.Hist.Surv_PB1(:,3) = fit(:,2);
            handles.Hist.Surv_PB1(:,4) = fit2(:,2);
            handles.Hist.Surv_PB1(:,5) = fit2(:,3);
            handles.Hist.Surv_PB1(:,6) = fit2(:,4);
            handles.FitPar.Surv_PB1 = [fitpar, fitpar2, PartialBoundFraction1; espSigma, espSigma2, BFerror1];
        end
    else
        TrLog1 = Hist1;
        fit = [Hist1(:,1),Hist1(:,3)];
        fit2 = [Hist1(:,1),Hist1(:,4:6)];
        if get(handles.PBCorrect,'Value') == 1
            fitpar = handles.FitPar.Surv_PB1(1,1:2);
            espSigma = handles.FitPar.Surv_PB1(2,1:2);
            fitpar2 = handles.FitPar.Surv_PB1(1,3:6);
            espSigma2 = handles.FitPar.Surv_PB1(2,3:6);
        else
            fitpar = handles.FitPar.Surv1(1,1:2);
            espSigma = handles.FitPar.Surv1(2,1:2);
            fitpar2 = handles.FitPar.Surv1(1,3:6);
            espSigma2 = handles.FitPar.Surv1(2,3:6);
        end
    end
    
%     HistAxes = findobj ('Tag', 'hist_axes');
    axes(handles.C1_C2SurvAxes);
    plot(TrLog1(:,1), TrLog1(:,2),'ok');
    hold on;
    plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,3),'b', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,4),'b--', 'LineWidth', 1);
    hold off
    xlabel('Time [s]');
    ylabel('Counts');
    title({['Fit of the survival probability',...
        ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction1,3), ' \pm ',...
        num2str(BFerror1,2)],...
        ['One Component fit k_{off} = ',num2str(fitpar(1),3),...
        ' \pm ', num2str(espSigma(1),2), ' s^{-1}; C_{eq} = ', num2str(fitpar(2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
        ' \pm ', num2str(espSigma2(1),2), ' s^{-1}'], ...
        ['k_{off, 2} = ', num2str(fitpar2(2),3),...
        ' \pm ', num2str(espSigma2(2),2), ' s^{-1}; f_1 = ', num2str(fitpar2(3),3),...
        '\pm', num2str(espSigma2(3),2), 'C_{eq} = ',  num2str(fitpar2(4),2)]});
    legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit');
    box on;
    
    
        
    
    
end
handles.CumHist1 = Hist1;
set(handles.FstatText1,'String',['F stat: ', num2str(Fstat1,4)]);
set(handles.pValText1,'String',['p-value: ', num2str(pval1,4)]);

%Channel 2
PartialBoundFraction2 = handles.PartBoundFrac2;
BFerror2 = handles.BFerror2;
if get(handles.What2Plot,'Value') == 1;
    

    if get(handles.PBCorrect,'Value') == 1 
        Hist2 = handles.Hist.Res_PB2;
        pval2 = handles.FitPar.Res2(1,8);
        Fstat2 = handles.FitPar.Res2(2,8);
    else
        Hist2 = handles.Hist.Res2;
        pval2 = handles.FitPar.Res_PB2(1,8);
        Fstat2 = handles.FitPar.Res_PB2(2,8);
    end
    
    if size(Hist2,2) < 3

        % Bin the histogram according to the selected settings;

        binFactor = handles.Parameters.bin + 1;
        nBins2 = floor(length(Hist2(:,1))/binFactor);
        
            TimeBins2 = reshape(Hist2(1:nBins2*binFactor,1), binFactor, nBins2);
            CountsBins2 = reshape(Hist2(1:nBins2*binFactor,2), binFactor, nBins2);
        

        if binFactor ~= 1
            resTimeHist_Binned2 = [];
            resTimeHist_Binned2(:,1) = mean(TimeBins2);
            resTimeHist_Binned2(:,2) = sum(CountsBins2);
        else
            resTimeHist_Binned2 = Hist2;
        end

        % Remove values < 0 from the histogram
        idx = find(resTimeHist_Binned2(:,2) <= 0);
        resTimeHist_Binned2(idx,:) = [];
        % Fit the binned residence time histogram with exponentials;
        [fitpar,espSigma, fit]= ExpDecay_fit_resTime(resTimeHist_Binned2, 1 );
        [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit_resTime(resTimeHist_Binned2, [1 0.1]);
        if get(handles.PBCorrect,'Value') == 1 
            handles.Hist.Res_PB2(:, 1:2) = resTimeHist_Binned2;
            handles.Hist.Res_PB2(:,3) = fit(:,2);
            handles.Hist.Res_PB2(:,4) = fit2(:,2);
            handles.Hist.Res_PB2(:,5) = fit2(:,3);
            handles.Hist.Res_PB2(:,6) = fit2(:,4);
            
            handles.FitPar.Res_PB2 = [fitpar, fitpar2, PartialBoundFraction2; espSigma, espSigma2, BFerror2];
        else
            handles.Hist.Res2(:, 1:2) = resTimeHist_Binned2;
            handles.Hist.Res2(:,3) = fit(:,2);
            handles.Hist.Res2(:,4) = fit2(:,2);
            handles.Hist.Res2(:,5) = fit2(:,3);
            handles.Hist.Res2(:,6) = fit2(:,4);
            
            handles.FitPar.Res2 = [fitpar, fitpar2, PartialBoundFraction2; espSigma, espSigma2, BFerror2];
        end
        
    else
        resTimeHist_Binned2 = Hist2;
        fit = [Hist2(:,1),Hist2(:,3)];
        fit2 = [Hist2(:,1),Hist2(:,4:6)];
        if get(handles.PBCorrect,'Value') == 1 
            fitpar = handles.FitPar.Res_PB2(1,1:2);
            espSigma = handles.FitPar.Res_PB2(2,1:2);
            fitpar2 = handles.FitPar.Res_PB2(1,3:6);
            espSigma2 = handles.FitPar.Res_PB2(2,3:6);
        else
            fitpar = handles.FitPar.Res2(1,1:2);
            espSigma = handles.FitPar.Res2(2,1:2);
            fitpar2 = handles.FitPar.Res2(1,3:6);
            espSigma2 = handles.FitPar.Res2(2,3:6);
        end
    end
        
    
    
    
%     HistAxes = findobj ('Tag', 'hist_axes');
    axes(handles.C2_C1SurvAxes);
    bar(resTimeHist_Binned2(:,1), resTimeHist_Binned2(:,2),'EdgeColor', [0 0 0], ...
        'FaceColor',[0 0 0], 'BarWidth', 1, 'LineWidth', 0.75);
    hold on;
    plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,3),'b','LineWidth',1);
    plot(fit2(:,1),fit2(:,4),'b--','LineWidth',1);
    hold off
    xlabel('Residence Time [s]')
    ylabel('Counts');
    title({['Fit of the residence time histogram',...
        ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction2,3), ' \pm ',...
        num2str(BFerror2,2)],...
        ['One Component fit k_{off} = ',num2str(fitpar(1),3),...
        ' \pm ', num2str(espSigma(1),2), ' s^{-1}; C_{eq} = ', num2str(fitpar(2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
        ' \pm ', num2str(espSigma2(1),2), ' s^{-1}'], ...
        ['k_{off, 2} = ', num2str(fitpar2(2),3),...
        ' \pm ', num2str(espSigma2(2),2), ' s^{-1}; f_1 = ', num2str(fitpar2(3),3),...
        '\pm', num2str(espSigma2(3),2), ' C_{eq} = ',  num2str(fitpar2(4),2)]});
    legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit','Component 1', 'Component 2');
    box on;
    
    
    
    
    


else
    if get(handles.PBCorrect,'Value') == 1 
        Hist2 = handles.Hist.Surv_PB2;
        pval2 = handles.FitPar.Surv_PB2(1,8);
        Fstat2 = handles.FitPar.Surv_PB2(2,8);
    else
        Hist2 = handles.Hist.Surv2;
        pval2 = handles.FitPar.Surv2(1,8);
        Fstat2 = handles.FitPar.Surv2(2,8);
    end
    
    if size(Hist2,2) < 3
    
        % Logarithmic sampling of the survival prob
        if handles.Parameters.bin ~= 0
            
            
            TimePoints2 = logspace(log10(min(Hist2(:,1))),...
                log10(max(Hist2(:,1))),handles.Parameters.bin);
            
            TrLog2 =[];
            
            TimePoints2 = round(TimePoints2*1000)/1000;
            for i = 1:length(TimePoints2)
                
                idx = find(Hist2(:,1) <= TimePoints2(i)+0.0001,1,'last');
                TrLog2(i,:) = Hist2(idx,:);
            end
            
        else
            TrLog2 = Hist2;
        end
        
        tpoint1 = TrLog2(1,1);
        TrLog2(:,1) = TrLog2(:,1) - tpoint1;
        
        
        
        
        % Fit exponentials to the Survival probability
        [fitpar,espSigma, fit]= ExpDecay_fit(TrLog2, 1/(tpoint1*10));
        [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit(TrLog2, [1/(tpoint1*10) 1/(tpoint1*100)]);
        fit(:,1) = fit(:,1) + tpoint1;
        fit2(:,1) = fit2(:,1) + tpoint1;
        TrLog2(:,1) = TrLog2(:,1) + tpoint1;
        if get(handles.PBCorrect,'Value') == 0
            handles.Hist.Surv2(:, 1:2) = TrLog2;
            handles.Hist.Surv2(:,3) = fit(:,2);
            handles.Hist.Surv2(:,4) = fit2(:,2);
            handles.Hist.Surv2(:,5) = fit2(:,3);
            handles.Hist.Surv2(:,6) = fit2(:,4);
            handles.FitPar.Surv2 = [fitpar, fitpar2, PartialBoundFraction2; espSigma, espSigma2, BFerror2];
        else
            handles.Hist.Surv_PB2(:, 1:2) = TrLog2;
            handles.Hist.Surv_PB2(:,3) = fit(:,2);
            handles.Hist.Surv_PB2(:,4) = fit2(:,2);
            handles.Hist.Surv_PB2(:,5) = fit2(:,3);
            handles.Hist.Surv_PB2(:,6) = fit2(:,4);
            handles.FitPar.Surv_PB2 = [fitpar, fitpar2, PartialBoundFraction2; espSigma, espSigma2, BFerror2];
        end
    else
        TrLog2 = Hist2;
        fit = [Hist2(:,1),Hist2(:,3)];
        fit2 = [Hist2(:,1),Hist2(:,4:6)];
        if get(handles.PBCorrect,'Value') == 1
            fitpar = handles.FitPar.Surv_PB2(1,1:2);
            espSigma = handles.FitPar.Surv_PB2(2,1:2);
            fitpar2 = handles.FitPar.Surv_PB2(1,3:6);
            espSigma2 = handles.FitPar.Surv_PB2(2,3:6);
        else
            fitpar = handles.FitPar.Surv2(1,1:2);
            espSigma = handles.FitPar.Surv2(2,1:2);
            fitpar2 = handles.FitPar.Surv2(1,3:6);
            espSigma2 = handles.FitPar.Surv2(2,3:6);
        end
    end
    
%     HistAxes = findobj ('Tag', 'hist_axes');
    axes(handles.C2_C1SurvAxes);
    plot(TrLog2(:,1), TrLog2(:,2),'ok');
    hold on;
    plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,3),'b', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,4),'b--', 'LineWidth', 1);
    hold off
    xlabel('Time [s]');
    ylabel('Counts');
    title({['Fit of the survival probability',...
        ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction2,3), ' \pm ',...
        num2str(BFerror2,2)],...
        ['One Component fit k_{off} = ',num2str(fitpar(1),3),...
        ' \pm ', num2str(espSigma(1),2), ' s^{-1}; C_{eq} = ', num2str(fitpar(2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
        ' \pm ', num2str(espSigma2(1),2), ' s^{-1}'], ...
        ['k_{off, 2} = ', num2str(fitpar2(2),3),...
        ' \pm ', num2str(espSigma2(2),2), ' s^{-1}; f_1 = ', num2str(fitpar2(3),3),...
        '\pm', num2str(espSigma2(3),2), 'C_{eq} = ',  num2str(fitpar2(4),2)]});
    legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit');
    box on;
    
    
        
    
    
end
handles.CumHist2 = Hist2;
set(handles.FstatText2,'String',['F stat: ', num2str(Fstat2,4)]);
set(handles.pValText2,'String',['p-value: ', num2str(pval2,4)]);
guidata(hObject,handles);

% --- Executes on button press in CopyHist.
function CopyHist_Callback(hObject, eventdata, handles)
% hObject    handle to CopyHist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.What2Plot,'Value') == 1 && ~get(handles.PBCorrect,'Value')
    Hist1 = handles.Hist.Res1;
    Hist2 = handles.Hist.Res2;
elseif get(handles.What2Plot,'Value') == 1 && get(handles.PBCorrect,'Value')
    Hist1 = handles.Hist.Res_PB1;
    Hist2 = handles.Hist.Res_PB2;
elseif get(handles.What2Plot,'Value') == 2 && ~get(handles.PBCorrect,'Value')
    Hist1 = handles.Hist.Surv1;
    Hist2 = handles.Hist.Surv2;
elseif get(handles.What2Plot,'Value') == 2 && get(handles.PBCorrect,'Value')
    Hist1 = handles.Hist.Surv_PB1;
    Hist2 = handles.Hist.Surv_PB2;
end
if size(Hist1,1) > size(Hist2,1)
    Hist2(size(Hist1,1),:) = zeros(1,size(Hist2,2));
elseif size(Hist2,1) > size(Hist1,1)
    Hist1(size(Hist2,1),:) = zeros(1,size(Hist1,2));
end
Hist = [Hist1, Hist2];
    
num2clip(Hist);

% --- Executes on button press in CopyAnalysisPar.
function CopyAnalysisPar_Callback(hObject, eventdata, handles)
% hObject    handle to CopyAnalysisPar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Prepare the vector output.
OutPar(1) = handles.Flag;
OutPar(2) = handles.Parameters.ThreshL;
OutPar(3) = handles.Parameters.ThreshH;
OutPar(4) = handles.Parameters.minBoundFrames;
OutPar(5) = handles.Parameters.bin;
OutPar(6) = handles.FrameTime;

num2clip(OutPar);


% --- Executes on button press in CopyFitPar.
function CopyFitPar_Callback(hObject, eventdata, handles)
% hObject    handle to CopyFitPar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.What2Plot,'Value') == 1 && ~get(handles.PBCorrect,'Value')
    FitPar1 = handles.FitPar.Res1;
    FitPar2 = handles.FitPar.Res2;
elseif get(handles.What2Plot,'Value') == 1 && get(handles.PBCorrect,'Value')
    FitPar1 = handles.FitPar.Res_PB1;
    FitPar2 = handles.FitPar.Res_PB2;
elseif get(handles.What2Plot,'Value') == 2 && ~get(handles.PBCorrect,'Value')
    FitPar1 = handles.FitPar.Surv1;
    FitPar2 = handles.FitPar.Surv2;
elseif get(handles.What2Plot,'Value') == 2 && get(handles.PBCorrect,'Value')
    FitPar1 = handles.FitPar.Surv_PB1;
    FitPar2 = handles.FitPar.Surv_PB2;
end

FitPar = [FitPar1; FitPar2];

num2clip(FitPar);

% --- Executes on button press in SaveFig.
function SaveFig_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Prepare data for figure;

Arrivals =  handles.ArrivalHist;
Departures =  handles.DepartHist;


if get(handles.What2Plot,'Value') == 1 && ~get(handles.PBCorrect,'Value')
    Hist1 = handles.Hist.Res1;
    FitParAll1 = handles.FitPar.Res1;
    Hist2 = handles.Hist.Res2;
    FitParAll2 = handles.FitPar.Res2;
elseif get(handles.What2Plot,'Value') == 1 && get(handles.PBCorrect,'Value')
    Hist1 = handles.Hist.Res_PB1;
    FitParAll1 = handles.FitPar.Res_PB1;
    Hist2 = handles.Hist.Res_PB2;
    FitParAll2 = handles.FitPar.Res_PB2;
elseif get(handles.What2Plot,'Value') == 2 && ~get(handles.PBCorrect,'Value')
    Hist1 = handles.Hist.Surv1;
    FitParAll1 = handles.FitPar.Surv1;
    Hist2 = handles.Hist.Surv2;
    FitParAll2 = handles.FitPar.Surv2;
elseif get(handles.What2Plot,'Value') == 2 && get(handles.PBCorrect,'Value')
    Hist1 = handles.Hist.Surv_PB1;
    FitParAll1 = handles.FitPar.Surv_PB1;
    Hist2 = handles.Hist.Surv_PB2;
    FitParAll2 = handles.FitPar.Surv_PB2;
end


% Hist = handles.Hist;
Flag = get(handles.What2Plot,'Value');

fitpar1 = FitParAll1(1, 1:2);
fitpar2 = FitParAll2(1, 1:2);
fitpar1_2 = FitParAll1(1, 3:6);
fitpar2_2 = FitParAll2(1, 3:6);
PartialBoundFraction1 = FitParAll1(1,7);
PartialBoundFraction2 = FitParAll2(1,7);

espSigma1 = FitParAll1(2, 1:2);
espSigma2 = FitParAll2(2, 1:2);
espSigma1_2 = FitParAll1(2, 3:6);
espSigma2_2 = FitParAll2(2, 3:6);
BFerror1 = FitParAll1(2, 7);
BFerror2 = FitParAll2(2, 7);




% Plot the photobleaching curve
scrsz = get(0,'ScreenSize');
h1 = figure('Position',[.15*scrsz(3) .15*scrsz(4) .5*scrsz(3) .4*scrsz(4)]);



subplot(2,2,1);
bar(Arrivals(:,1),Arrivals(:,2));
title('Relative Arrival');
xlabel('Time (s)');
ylabel('Frequency');

subplot(2,2,2);
bar(Departures(:,1),Departures(:,2));
title('Relative Departures');
xlabel('Time (s)');
ylabel('Frequency');
box on;

% Plot the residence time/ cumulative histogram.

subplot(2,2,3);
plot(Hist1(:,1), Hist1(:,2),'ok');
hold on;
plot(Hist1(:,1),Hist1(:,3),'g', 'LineWidth', 1);
plot(Hist1(:,1),Hist1(:,4),'r', 'LineWidth', 1);
plot(Hist1(:,1),Hist1(:,5),'b', 'LineWidth', 1);
plot(Hist1(:,1),Hist1(:,6),'b--', 'LineWidth', 1);
hold off;
box on;
if Flag == 1
    xlabel('Residence Time [s]')
    ylabel('Counts');
    title({['Fit of the residence time histogram',...
        ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction1,3), ' \pm ',...
        num2str(BFerror1,2)],...
        ['One Component fit k_{off} = ',num2str(fitpar1(1),3),...
        ' \pm ', num2str(espSigma1(1),2), ' s; C_eq = ', num2str(fitpar1(2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(fitpar1_2(1),3),...
        ' \pm ', num2str(espSigma1_2(1),2), 's'], ...
        ['k_{off, 2} = ', num2str(fitpar1_2(2),3),...
        ' \pm ', num2str(espSigma1_2(2),2), 's; f_1 = ', num2str(fitpar1_2(3),3),...
        '\pm', num2str(espSigma1_2(3),2), 'C_eq = ',  num2str(fitpar1_2(4),2)]});
    legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit','Component 1','Component 2');
    box on;
    DefaultName = handles.PathName;
    DefaultName = [DefaultName, 'GlobalResTHist'];
    
elseif Flag == 2;
    xlabel('Time [s]');
    ylabel('Counts');
    title({['Fit of the survival probability',...
        ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction1,3), ' \pm ',...
        num2str(BFerror1,2)],...
        ['One Component fit k_{off} = ',num2str(fitpar1(1),3),...
        ' \pm ', num2str(espSigma1(1),2), ' s; C_eq = ', num2str(fitpar1(2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(fitpar1_2(1),3),...
        ' \pm ', num2str(espSigma1_2(1),2), 's'], ...
        ['k_{off, 2} = ', num2str(fitpar1_2(2),3),...
        ' \pm ', num2str(espSigma1_2(2),2), 's; f_1 = ', num2str(fitpar1_2(3),3),...
        '\pm', num2str(espSigma1_2(3),2), 'C_eq = ',  num2str(fitpar1_2(4),2)]});
    legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit');
    box on;
    DefaultName = handles.PathName;
    DefaultName = [DefaultName, 'GlobalSurvTHist'];
    
    
end
subplot(2,2,4);
plot(Hist2(:,1), Hist2(:,2),'ok');
hold on;
plot(Hist2(:,1),Hist2(:,3),'g', 'LineWidth', 1);
plot(Hist2(:,1),Hist2(:,4),'r', 'LineWidth', 1);
plot(Hist2(:,1),Hist2(:,5),'b', 'LineWidth', 1);
plot(Hist2(:,1),Hist2(:,6),'b--', 'LineWidth', 1);
hold off;
box on;
if Flag == 1
    xlabel('Residence Time [s]')
    ylabel('Counts');
    title({['Fit of the residence time histogram',...
        ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction2,3), ' \pm ',...
        num2str(BFerror2,2)],...
        ['One Component fit k_{off} = ',num2str(fitpar2(1),3),...
        ' \pm ', num2str(espSigma2(1),2), ' s; C_eq = ', num2str(fitpar2(2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(fitpar2_2(1),3),...
        ' \pm ', num2str(espSigma2_2(1),2), 's'], ...
        ['k_{off, 2} = ', num2str(fitpar2_2(2),3),...
        ' \pm ', num2str(espSigma2_2(2),2), 's; f_1 = ', num2str(fitpar2_2(3),3),...
        '\pm', num2str(espSigma2_2(3),2), 'C_eq = ',  num2str(fitpar2_2(4),2)]});
    legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit','Component 1','Component 2');
    box on;
    DefaultName = handles.PathName;
    DefaultName = [DefaultName, 'GlobalResTHist'];
    
elseif Flag == 2;
    xlabel('Time [s]');
    ylabel('Counts');
    title({['Fit of the survival probability',...
        ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction2,3), ' \pm ',...
        num2str(BFerror2,2)],...
        ['One Component fit k_{off} = ',num2str(fitpar2(1),3),...
        ' \pm ', num2str(espSigma2(1),2), ' s; C_eq = ', num2str(fitpar2(2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(fitpar2_2(1),3),...
        ' \pm ', num2str(espSigma2_2(1),2), 's'], ...
        ['k_{off, 2} = ', num2str(fitpar2_2(2),3),...
        ' \pm ', num2str(espSigma2_2(2),2), 's; f_1 = ', num2str(fitpar2_2(3),3),...
        '\pm', num2str(espSigma2_2(3),2), 'C_eq = ',  num2str(fitpar2_2(4),2)]});
    legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit');
    box on;
    DefaultName = handles.PathName;
    DefaultName = [DefaultName, 'GlobalSurvTHist'];
    
    
end
%Save Figure

FilterSpec = {'*.fig';'*.ai';'*.bpm';'*.eps';'*.jpeg';'*.pdf';'*.tif'};
[FileNameOut,PathNameOut] = uiputfile(FilterSpec,'Save Figure',DefaultName);
if FileNameOut ~= 0;
    saveas(h1, [PathNameOut,FileNameOut]);
end


% --- Executes on button press in BoxPlotButton.
function BoxPlotButton_Callback(hObject, eventdata, handles)
% hObject    handle to BoxPlotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in CopyArrDepart.
function CopyArrDepart_Callback(hObject, eventdata, handles)
% hObject    handle to CopyArrDepart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Arrivals =  handles.ArrivalHist;
Departures =  handles.DepartHist;

num2clip([Arrivals, Departures]);


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SaveData_Callback(hObject, eventdata, handles)
% hObject    handle to SaveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'PathName')
    defPath = handles.PathName;
else
    defPath = pwd;
end

[fname,pname] = uiputfile('*.mat','Save Residence Time Analysis',[defPath filesep '2color_Binding.mat']);

if fname ~= 0
    Results.Hist = handles.Hist;
    Results.FitPar = handles.FitPar;
    Results.NParticles1 = handles.NParticles1;
    Results.BleachRates1 = handles.BleachRates1;
    Results.NParticles2 = handles.NParticles2;
    Results.BleachRates2 = handles.BleachRates2;
    Results.Parameters = handles.Parameters;
    Results.totalTracks1 = handles.totalTracks1;
    Results.totalTracks2 = handles.totalTracks2;
    Results.BoundMolecules1 = handles.BoundMolecules1;
    Results.TotalMolecules1 = handles.TotalMolecules1;
    Results.BoundMolecules2 = handles.BoundMolecules2;
    Results.TotalMolecules2 = handles.TotalMolecules2;
    Results.totalROIs = handles.nROIs;
    Results.ArrivalHist =  handles.ArrivalHist;
    Results.DepartHist =  handles.DepartHist;
    Results.ColocParticles = handles.ColocParticles;
    Results.coboundTracks1_all = handles.coboundTracks1_all;
    Results.coboundTracks2_all = handles.coboundTracks2_all;
    Results.ArrivalTimes_rel = handles.ArrivalTimes_rel;
    Results.DepartTimes_rel = handles.DepartTimes_rel;
    Results.TrkNum_ArrDep = handles.TrkNum_ArrDep;
    save([pname,fname],'Results');
end


% --------------------------------------------------------------------
function LoadData_Callback(hObject, eventdata, handles)
% hObject    handle to LoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'PathName')
    defPath = handles.PathName;
else
    defPath = pwd;
end
[fname,pname] = uigetfile('*.mat','Load Residence Time Analysis',defPath );

if fname ~= 0
    IN = load(fullfile(pname,fname));
    handles.PathName = pname;
    handles.Hist = IN.Results.Hist;
    handles.FitPar = IN.Results.FitPar;
    handles.NParticles1 = IN.Results.NParticles1;
    handles.NParticles2 = IN.Results.NParticles2;
    handles.BleachRates1 = IN.Results.BleachRates1;
    handles.BleachRates2 = IN.Results.BleachRates2;
    handles.Parameters = IN.Results.Parameters;
    handles.totalTracks1 = IN.Results.totalTracks1;
    handles.totalTracks2 = IN.Results.totalTracks2;
    handles.BoundMolecules1 = IN.Results.BoundMolecules1;
    handles.BoundMolecules2 = IN.Results.BoundMolecules2;
    handles.TotalMolecules1 = IN.Results.TotalMolecules1;
    handles.TotalMolecules2 = IN.Results.TotalMolecules2;
    handles.CumHist1 = handles.Hist.Surv_PB1;
    handles.CumHist2 = handles.Hist.Surv_PB2;
    handles.ArrivalHist = IN.Results.ArrivalHist;
    handles.DepartHist = IN.Results.DepartHist;
    handles.ColocParticles = IN.Results.ColocParticles;
    if isfield(IN.Results,'totalROIs');
        handles.nROIs = IN.Results.totalROIs;
        set(handles.totalROIsText,'String',['Total ROIs: ',num2str(handles.nROIs)]);
    else
        set(handles.totalROIsText,'String','Total ROIs: Not Available');
    end
    ArrDep_title_flag = 0;
    if isfield(IN.Results, 'ArrivalTimes_rel')
        handles.ArrivalTimes_rel_all = IN.Results.ArrivalTimes_rel;
        handles.DepartTimes_rel_all = IN.Results.DepartTimes_rel;
        handles.TrkNum_ArrDep = IN.Results.TrkNum_ArrDep;
        ArrDep_title_flag = 1;
    end
    
    %Display Survival plot with PB correction by default
    set(handles.What2Plot,'Value',2);
    set(handles.PBCorrect,'Value',1);
    %Display track & particle counts
     set(handles.totalTracksText,'String',['Total Tracks: Ch1: ',num2str(handles.totalTracks1), ', Ch2: ',num2str(handles.totalTracks2) ]);
    set(handles.boundParticlesText,'String',['Bound Particles: Ch1: ',num2str(handles.BoundMolecules1), ', Ch2: ',num2str(handles.BoundMolecules2)]);
    set(handles.totalParticlesText,'String',['Total Particles: Ch1: ',num2str(handles.TotalMolecules1), ' Ch2: ',num2str(handles.TotalMolecules2)]);
    set(handles.totalColocPartText,'String',['Total Coloc Particles: ', num2str(handles.ColocParticles)]);
    
    %Display F stat info
    set(handles.FstatText1,'String',['F stat: ', num2str(handles.FitPar.Surv_PB1(2,8))]);
    set(handles.pValText1,'String', ['p-value: ', num2str(handles.FitPar.Surv_PB1(1,8))]);
    
    set(handles.FstatText2,'String',['F stat: ', num2str(handles.FitPar.Surv_PB2(2,8))]);
    set(handles.pValText2,'String', ['p-value: ', num2str(handles.FitPar.Surv_PB2(1,8))]);
    
    %Plot the arrivals and departures
    axes(handles.ArrivalAxes)
    bar(handles.ArrivalHist(:,1),handles.ArrivalHist(:,2));
    
    xlabel('Time (s)');
    ylabel('Frequency');
    box on;
    
    if ArrDep_title_flag == 1
        Arrive_mean = mean(handles.ArrivalTimes_rel_all);
        Arrive_median = median(handles.ArrivalTimes_rel_all);
        Arrive_std = std(handles.ArrivalTimes_rel_all);
        Arrive_mode = mode(handles.ArrivalTimes_rel_all);
        
        title({['Relative Arrival - mean = ', num2str(Arrive_mean,2), ', std. dev. = ', num2str(Arrive_std,2)],...
        ['median = ', num2str(Arrive_median,2), ', mode = ',num2str(Arrive_mode,2), ', N = ',num2str(size(handles.ArrivalTimes_rel_all,1))]});
    else
        title('Relative Arrival');
    end
    
    axes(handles.DepartAxes);
    bar(handles.DepartHist(:,1),handles.DepartHist(:,2));
    if ArrDep_title_flag == 1
        Depart_mean = mean(handles.DepartTimes_rel_all);
        Depart_median = median(handles.DepartTimes_rel_all);
        Depart_std = std(handles.DepartTimes_rel_all);
        Depart_mode = mode(handles.DepartTimes_rel_all);
        
        title({['Relative Arrival - mean = ', num2str(Depart_mean,2), ', std. dev. = ', num2str(Depart_std,2)],...
        ['median = ', num2str(Depart_median,2), ', mode = ',num2str(Depart_mode,2), ', N = ',num2str(size(handles.DepartTimes_rel_all,1))]});
    else
        title('Relative Departure');
    end
    xlabel('Time (s)');
    ylabel('Frequency');
    box on;
    
   
    
    %Plot the survival distribution and fits
    %Channel 1
    axes(handles.C1_C2SurvAxes)
    plot(handles.Hist.Surv_PB1(:,1), handles.Hist.Surv_PB1(:,2),'ok');
    hold on;
    plot(handles.Hist.Surv_PB1(:,1),handles.Hist.Surv_PB1(:,3),'g', 'LineWidth', 1);
    plot(handles.Hist.Surv_PB1(:,1),handles.Hist.Surv_PB1(:,4),'r', 'LineWidth', 1);
    plot(handles.Hist.Surv_PB1(:,1),handles.Hist.Surv_PB1(:,5),'b', 'LineWidth', 1);
    plot(handles.Hist.Surv_PB1(:,1),handles.Hist.Surv_PB1(:,6),'b--', 'LineWidth', 1);
    hold off;
    xlabel('Time [s]');
    ylabel('Counts');
    title({['Fit of the survival probability',...
        ' - Measured [partial] C_{eq} = ' num2str(handles.FitPar.Surv_PB1(1,7),3), ' \pm ',...
        num2str(handles.FitPar.Surv_PB1(2,7),2)],...
        ['One Component fit k_{off} = ',num2str(handles.FitPar.Surv_PB1(1,1),3),...
        ' \pm ', num2str(handles.FitPar.Surv_PB1(2,7),2), ' s^{-1}; C_{eq} = ', num2str(handles.FitPar.Surv_PB1(1,2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(handles.FitPar.Surv_PB1(1,3),3),...
        ' \pm ', num2str(handles.FitPar.Surv_PB1(2,3),2), ' s^{-1}'], ...
        ['k_{off, 2} = ', num2str(handles.FitPar.Surv_PB1(1,4),3),...
        ' \pm ', num2str(handles.FitPar.Surv_PB1(2,4),2), ' s^{-1}; f_1 = ', num2str(handles.FitPar.Surv_PB1(1,5),3),...
        '\pm', num2str(handles.FitPar.Surv_PB1(2,5),2), 'C_{eq} = ',  num2str(handles.FitPar.Surv_PB1(1,6),2)]});
    legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit');
    box on;
    
    %Channel 2
    axes(handles.C2_C1SurvAxes)
    plot(handles.Hist.Surv_PB2(:,1), handles.Hist.Surv_PB2(:,2),'ok');
    hold on;
    plot(handles.Hist.Surv_PB2(:,1),handles.Hist.Surv_PB2(:,3),'g', 'LineWidth', 1);
    plot(handles.Hist.Surv_PB2(:,1),handles.Hist.Surv_PB2(:,4),'r', 'LineWidth', 1);
    plot(handles.Hist.Surv_PB2(:,1),handles.Hist.Surv_PB2(:,5),'b', 'LineWidth', 1);
    plot(handles.Hist.Surv_PB2(:,1),handles.Hist.Surv_PB2(:,6),'b--', 'LineWidth', 1);
    hold off;
    xlabel('Time [s]');
    ylabel('Counts');
    title({['Fit of the survival probability',...
        ' - Measured [partial] C_{eq} = ' num2str(handles.FitPar.Surv_PB2(1,7),3), ' \pm ',...
        num2str(handles.FitPar.Surv_PB2(2,7),2)],...
        ['One Component fit k_{off} = ',num2str(handles.FitPar.Surv_PB2(1,1),3),...
        ' \pm ', num2str(handles.FitPar.Surv_PB2(2,7),2), ' s^{-1}; C_{eq} = ', num2str(handles.FitPar.Surv_PB2(1,2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(handles.FitPar.Surv_PB2(1,3),3),...
        ' \pm ', num2str(handles.FitPar.Surv_PB2(2,3),2), ' s^{-1}'], ...
        ['k_{off, 2} = ', num2str(handles.FitPar.Surv_PB2(1,4),3),...
        ' \pm ', num2str(handles.FitPar.Surv_PB2(2,4),2), ' s^{-1}; f_1 = ', num2str(handles.FitPar.Surv_PB2(1,5),3),...
        '\pm', num2str(handles.FitPar.Surv_PB2(2,5),2), 'C_{eq} = ',  num2str(handles.FitPar.Surv_PB2(1,6),2)]});
    legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit');
    box on;
    
    
    %Display Parameters
    set(handles.frame2frameMaxText,'String',['Max. Frame to Frame Jump: ', num2str(handles.Parameters.ThreshL), ' um']);
    set(handles.End2EndMaxText,'String',['Max. End to End Jump: ', num2str(handles.Parameters.ThreshH), ' um']);
    set(handles.NminText,'String',['Min. Bound Frames (Nmin): ', num2str(handles.Parameters.minBoundFrames)]);
    guidata(hObject,handles);
end
