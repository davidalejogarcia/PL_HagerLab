function varargout = AnalyzeResTimes_GUI(varargin)
% ANALYZERESTIMES_GUI MATLAB code for AnalyzeResTimes_GUI.fig
%      ANALYZERESTIMES_GUI, by itself, creates a new ANALYZERESTIMES_GUI or raises the existing
%      singleton*.
%
%      H = ANALYZERESTIMES_GUI returns the handle to a new ANALYZERESTIMES_GUI or the handle to
%      the existing singleton*.
%
%      ANALYZERESTIMES_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYZERESTIMES_GUI.M with the given input arguments.
%
%      ANALYZERESTIMES_GUI('Property','Value',...) creates a new ANALYZERESTIMES_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AnalyzeResTimes_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AnalyzeResTimes_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AnalyzeResTimes_GUI

% Last Modified by GUIDE v2.5 22-Sep-2017 11:33:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @AnalyzeResTimes_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @AnalyzeResTimes_GUI_OutputFcn, ...
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


% --- Executes just before AnalyzeResTimes_GUI is made visible.
function AnalyzeResTimes_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AnalyzeResTimes_GUI (see VARARGIN)

% Choose default command line output for AnalyzeResTimes_GUI
handles.output = hObject;
tic;
% Import Parameters.
if ~isempty(varargin)
    handles.Flag = varargin{1};
    handles.PathName = varargin{2};
    handles.FileNames = varargin{3};
    AnalysisS = varargin{4};

    handles.Parameters = [];
    handles.Parameters.ThreshL = str2num(AnalysisS{1});
    handles.Parameters.ThreshH = str2num(AnalysisS{2});
    handles.Parameters.minBoundFrames = str2num(AnalysisS{3});
    handles.Parameters.bin = str2num(AnalysisS{4});
    handles.FrameTime = str2num(AnalysisS{5});
    handles.maxTime = str2num(AnalysisS{6});




    % If only one file has been selected change the fileNames variable
    % to a cell array

    if ~iscell(handles.FileNames)
        handles.FileNames = {handles.FileNames};
    end;
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

        if isfield(Temp.Results, 'PreAnalysis') && ...
                isfield(Temp.Results.PreAnalysis, 'Tracks_um');
            if handles.maxTime > 0
                Tracks{i} = Temp.Results.PreAnalysis.Tracks_um(Temp.Results.PreAnalysis.Tracks_um(:,3) <= handles.maxTime,:);
                Particles{i} = Temp.Results.Tracking.Particles(Temp.Results.Tracking.Particles(:,6) <= handles.maxTime,:);
                NParticles{i} = Temp.Results.PreAnalysis.NParticles(Temp.Results.PreAnalysis.NParticles(:,1) <= handles.maxTime,:);
            else
                Tracks{i} = Temp.Results.PreAnalysis.Tracks_um;
                Particles{i} = Temp.Results.Tracking.Particles;
                NParticles{i} = Temp.Results.PreAnalysis.NParticles;
            end
            if isfield(Temp.Results,'Process') %check that the strcutures exist
                if isfield(Temp.Results.Process,'ROIlabel')
                    if isfield(Temp.Results.Process,'AllROIClasses') && ROIorClass == -1
                        ROIClass =  questdlg('Do you want to separate data based on ROI names or Class Names?','ROI or Class','ROI','Class','Class');
                        if strcmp(ROIClass,'ROI');
                            ROIorClass = 0;
                        else
                            ROIorClass = 1;
                        end
                    elseif ~isfield(Temp.Results.Process,'AllROIClasses')
                        ROIorClass = 0;
                    end
                    if ROIorClass == 0
                        roiLabels = Temp.Results.Process.ROIlabel;
                    elseif ROIorClass == 1
                        if ~isfield(Temp.Results.Process,'AllROIClasses')
                            errordlg(['File does not contain Class data: ,' handles.FileNames{i}]);
                        else
                            roiLabels = Temp.Results.Process.AllROIClasses;
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
                                if size(Tracks{i},2) < 8
                                    Tracks{i} = Tracks{i}(Tracks{i}(:,5) == ROIidx,:);
                                elseif size(Tracks{i},2) == 8
                                    Tracks{i} = Tracks{i}(Tracks{i}(:,6) == ROIidx,:);
                                else
                                    Tracks{i} = Tracks{i}(Tracks{i}(:,9) == ROIidx,:);
                                end
                                NParticles{i} = [NParticles{i}(:,1) NParticles{i}(:,ROIidx+1)];
                                ROIChoice{end+1,1} = Temp.Results.Process.ROIlabel{ROIidx};
                                nROIs = nROIs + 1;
                                
                            else
                                roiAct = [];
                                tracks_tmp = Tracks{i};
                                Tracks{i} = [];
                                for m = 1:size(Temp.Results.Process.ROIClass);
                                    if strcmpi(roiLabels{ROIidx,:},Temp.Results.Process.ROIClass{m,:})
                                        if size(Tracks{i},2) < 8
                                            Tracks{i} = [Tracks{i}; tracks_tmp(tracks_tmp(:,5) == m,:)];
                                        elseif size(Tracks{i},2) == 8
                                            Tracks{i} = [Tracks{i}; tracks_tmp(tracks_tmp(:,6) == m,:)];
                                        else
                                            Tracks{i} = [Tracks{i}; tracks_tmp(tracks_tmp(:,9) == m,:)];
                                        end
                                        roiAct = [roiAct;m];
                                    end

                                end
                                tmp = NParticles{i}(:,roiAct+1);
                                nROIs = nROIs + size(tmp,2);

                                NParticles{i} = [NParticles{i}(:,1) sum(tmp,2)];
                                ROIChoice{end+1,1} = Temp.Results.Process.AllROIClasses{ROIidx};
                            end

                        else
                            tmp = NParticles{i}(:,2:size(Temp.Results.Process.ROIlabel,1)+1);
                            nROIs = nROIs + size(tmp,2);
                            NParticles{i} = [NParticles{i}(:,1) sum(tmp,2)];
                            ROIChoice{end+1,1} = 'All';
                        end
                    else
                        ROIChoice{1} = roiLabels{1};
                        if ROIorClass == 1
                            tmp = NParticles{i}(:,2:size(NParticles{i},2));
                            NParticles{i} = [NParticles{i}(:,1) sum(tmp,2)];
                            nROIs = nROIs + size(tmp,2);
                        end
                    end
                else
                    ROIChoice{1} = '';
                end
            else
                ROIChoice{1} = '';
            end


            if ~isempty(Tracks{i})
                if size(Tracks{i},2) <= 8
                    [ImmTracks, Dummy] = calculateImmobileTracks...
                        (Tracks{i}, handles.Parameters.ThreshL,...
                        handles.Parameters.minBoundFrames, handles.Parameters.ThreshH, 0);
    

                    if isempty(ImmTracks)
                        TrackLengthHist{i} = [0 0];
                    else

                        % Identify bound molecules for each of the movies;
%                         [DiffCoef{i}, ResTimes{i}] = diffusionVersusResTime(ImmTracks,handles.FrameTime);
                        % Calculate survival histogram;
                        TrackLengthHist{i} = ...
                            calculateTrackLength(ImmTracks, ...
                            handles.FrameTime,handles.Parameters.minBoundFrames);



                    end
                else
                    ImmTracks = Tracks{i}(Tracks{i}(:,2) >= handles.Parameters.minBoundFrames,:);
                    ImmTracks = ImmTracks(ImmTracks(:,3) <= handles.Parameters.ThreshL,:);
                    if ~isempty(ImmTracks)
                        TrackLength = ImmTracks(:,2);
                        LongestTrack = max(TrackLength);
                        ShortestTrack = handles.Parameters.minBoundFrames;

                        TrackLengthHist{i} = zeros(LongestTrack-ShortestTrack + 1,2);
                        TrackLengthHist{i}(:,1) = (ShortestTrack-1: LongestTrack)* handles.FrameTime;
                        %         TrackLengthHist(:,2) = hist(TrackLength,(ShortestTrack:LongestTrack));
                        for  j = 1:length(TrackLength);
                            TrackLengthHist{i}(1:TrackLength(j)- ShortestTrack + 1,2) = ...
                                TrackLengthHist{i}(1:TrackLength(j) - ShortestTrack + 1,2)+1;
                        end;
                    end
                end
            end


        else
            errordlg(['Some of the files have not been preprocessed: ',handles.PathName,handles.FileNames{i}]);
            return
        end
    end
    handles.nROIs = nROIs;
    %Calculate the avg image intensity over time
    % avgInt = sum(pxIntesity)./sum(npx);

    % Calculate the total number of particles;

    % find the data containing the lowest number of frames;
    TimePoints = [];
    
%     DiffCoef_all = DiffCoef{1};
%     ResTimes_all = ResTimes{1};
%     for i = 2:nFiles
%         DiffCoef_all = [DiffCoef_all; DiffCoef{i}];
%         ResTimes_all = [ResTimes_all; ResTimes{i}];
%     end
    

    nMax = length(NParticles{1}(:,1));
    TimePoints = NParticles{1}(:,1);
    for i = 1:nFiles
        if isempty(NParticles{i})
            blag = 0;
        end
        if size(NParticles{i},1) < nMax;

            TimePoints = NParticles{i}(:,1);
            nMax = length(TimePoints);
        end
    end

    if isempty(TimePoints)
        errordlg('Too few bound particles to produce an histogram');
        return
    end

    % Accumulate the histogram
    Hist_Matrix = zeros(nMax,nFiles);

    for i = 1:nFiles
        Hist_Matrix(1:nMax,i) = NParticles{i}(1:nMax,2);
    end
    CumNParticles(:,1) = (TimePoints-1) * handles.FrameTime;
    CumNParticles(:,2) = sum(Hist_Matrix, 2);

    %Calculate the avg intensity


    % Accumulate the different survival histograms;

    % find the histogram containing the longest track
    nMax = 1;
    TimePoints = [];
    totalTracks = 0;
    for i = 1:nFiles
        trInd = unique(Tracks{i}(:,4));
        totalTracks = totalTracks + length(trInd);
        if isempty(TrackLengthHist{i})
            blag = 0;
        end
        if size(TrackLengthHist{i},1) > nMax;

            TimePoints = TrackLengthHist{i}(:,1);
            nMax = length(TimePoints);
        end
    end
    handles.totalTracks = totalTracks;
    if isempty(TimePoints)
        errordlg('Too few bound particles to produce an histogram');
        return
    end


    % Accumulate the histogram
    Hist_Matrix = zeros(nMax,nFiles);

    for i = 1:nFiles
        if ~isempty(TrackLengthHist{i})
            lengthHist = length(TrackLengthHist{i}(:,1));
            Hist_Matrix(1:lengthHist,i) = TrackLengthHist{i}(:,2);
        end
    end
    save('Survival_IndividualExperiments','Hist_Matrix','TimePoints')

    CumTrackLengthHist(:,1) = TimePoints;
    CumTrackLengthHist(:,2) = sum(Hist_Matrix, 2);
    

    % Copy Histogram to handles
    handles.CumHist = CumTrackLengthHist;


    % Calculate the fraction of bound molecules.

    % First calculate the residence time histogram
    ResTimeHist(:,1)= CumTrackLengthHist(:,1);
    ResTimeHist(:,2) = ...
        [(CumTrackLengthHist(1:end-1,2) - CumTrackLengthHist(2:end,2));CumTrackLengthHist(end,2)];

    % Then calculate the total number of bound spots
    TotalBoundMolecules = sum(ResTimeHist(:,2).*ResTimeHist(:,1)/handles.FrameTime);
    handles.BoundMolecules = TotalBoundMolecules;
    % and the total number of spots
    % TotalMolecules = sum(CumNParticles(:,2));
    TotalMolecules = 0;
    for i = 1:nFiles
        TotalMolecules = TotalMolecules + sum(NParticles{i}(:,2));
    end
    handles.TotalMolecules = TotalMolecules;
    % Finally divide the two.
    PartialBoundFraction = TotalBoundMolecules/TotalMolecules;

    BFerror = (sqrt(TotalBoundMolecules)/TotalBoundMolecules + ...
        sqrt(TotalMolecules)/TotalMolecules)*PartialBoundFraction;

    %%%Edited David A. Garcia
    save('DwellTimeDistribution','ResTimeHist','CumNParticles')

    % NOTE: THIS IS ONLY THE PARTIAL BOUND FRACTION
    % BECAUSE WE ARE LOOKING ONLY AT PARTICLES BOUND
    % FOR MORE THAN Nmin. We HAVE TO FIT THE HISTOGRAM OF RESIDENCE
    % TIMES AND EXTRAPOLATE THE TRUE VALUE OF THE BOUND FRACTION

    % Photobleaching correction




    % Fit a double exponential to the photobleaching curve
    disp(' ')
    disp('____________________________________')
    disp('Estimating bleaching characteristics')
    disp('____________________________________')

    [BleachRates, Dummy, CumNParticles(:,3:6)] =...
        ExpDecay_2Cmp_fit(CumNParticles, [1 0.1]);
    disp(['Bleach Rate 1: ', num2str(BleachRates(1), 3), ' s^-1'])
    disp(['Bleach Rate 2: ', num2str(BleachRates(2), 3), ' s^-1'])
    disp(['Fraction 1: ', num2str(BleachRates(3), 3)])
    disp('____________________________________')
    disp(' ')

    CumNParticles = CumNParticles(:,1:4);
    % Plot the photobleaching curve
    PhotobAxes = findobj ('Tag', 'photob_axes');
    axes(PhotobAxes);






    plot(CumNParticles(:,1), CumNParticles(:,2)/max(CumNParticles(:,4)), '+', 'MarkerEdgeColor', [0.5 0.5 0.5]');
    hold on
    plot(CumTrackLengthHist(:,1), CumTrackLengthHist(:,2)/max(CumTrackLengthHist(:,2)),'ok');
    plot(CumNParticles(:,3), CumNParticles(:,4)/max(CumNParticles(:,4)),'r');
    hold off;

    box on;

    title({'Photobleaching kinetics:', ['k_{b1} = ', num2str(BleachRates(1), 3),...
        ' s^{-1}, k_{b2} = ', num2str(BleachRates(2), 3), ' s^{-1}, f_1 = ',  ...
        num2str(BleachRates(3), 3)]})

    hold off;
    legend('Photobleaching Decay', 'Bound molecules Decay');

    xlabel('Time [s]');
    ylabel('Normalized Counts');
    % Copy CumNparticles to handles
    handles.NParticles = CumNParticles;
    handles.BleachRates = BleachRates;

    CumTrackLengthHist_PB = CumTrackLengthHist;
    CumTrackLengthHist_PB(:,2) = CumTrackLengthHist_PB(:,2)./ ...
        (BleachRates(3)*exp(-BleachRates(1).* CumTrackLengthHist_PB(:,1)) + ...
        (1-BleachRates(3))*exp(-BleachRates(2).* CumTrackLengthHist_PB(:,1)));
    ResTimeHist_PB(:,1)= CumTrackLengthHist_PB(:,1);
    ResTimeHist_PB(:,2) = ...
        [(CumTrackLengthHist_PB(1:end-1,2) - CumTrackLengthHist_PB(2:end,2));CumTrackLengthHist_PB(end,2)];

    %only use values where there are particles with that residence time for curve-fitting
    ind2fit = find(ResTimeHist(:,2) > 0);


    % Ask if you want to correct the histogram for photobleaching;
%     answer = questdlg('Do you want to correct the histogram for photobleaching?');
    answer = 'Yes';
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
    deltaT = ResTimeHist(2,1) - ResTimeHist(1,1);
    % Normalize the Survival probability for the partial bound fraction
    % if size(Tracks{1},2) < 8
    CumTrackLengthHist(:,2) = CumTrackLengthHist(:,2)/CumTrackLengthHist(1,2)...
        .*PartialBoundFraction;

    CumTrackLengthHist_PB(:,2) = CumTrackLengthHist_PB(:,2)/CumTrackLengthHist_PB(1,2)...
        .*PartialBoundFraction;
    % Normalize the residence time histogram for the partial bound fraction
    ResTimeHist(:,2) = ResTimeHist(:,2)/(deltaT*sum(ResTimeHist(:,2)))*...
        PartialBoundFraction;

    ResTimeHist_PB(:,2) = ResTimeHist_PB(:,2)/(deltaT*sum(ResTimeHist_PB(:,2)))*...
        PartialBoundFraction;
    % end


    % handles.Hist.Surv(:,1:2) = CumTrackLengthHist;
    % 
    % handles.Hist.Surv_PB(:,1:2) = CumTrackLengthHist_PB;
    % 
    % handles.Hist.Res(:,1:2) = ResTimeHist;
    % 
    % handles.Hist.Res_PB(:,1:2) = ResTimeHist_PB;
    handles.PartBoundFrac = PartialBoundFraction;
    handles.BFerror = BFerror;

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
            binFactor = floor(length(ResTimeHist(:,1))/nBins);
    %     binFactor = handles.Parameters.bin + 1;
    %     nBins = floor(length(ResTimeHist(:,1))/binFactor);
    %     if strcmp(answer,'Yes')
            TimeBins = reshape(ResTimeHist(1:nBins*binFactor,1), binFactor, nBins);
            CountsBins = reshape(ResTimeHist(1:nBins*binFactor,2), binFactor, nBins);
    %     else
            TimeBins_PB = reshape(ResTimeHist_PB(1:nBins*binFactor,1), binFactor, nBins);
            CountsBins_PB = reshape(ResTimeHist_PB(1:nBins*binFactor,2), binFactor, nBins);
    %     end


            resTimeHist_Binned = [];
            resTimeHist_Binned(:,1) = mean(TimeBins);
            resTimeHist_Binned(:,2) = sum(CountsBins);
            resTimeHist_Binned_PB = [];
            resTimeHist_Binned_PB(:,1) = mean(TimeBins_PB);
            resTimeHist_Binned_PB(:,2) = sum(CountsBins_PB);
        else
            resTimeHist_Binned = ResTimeHist;
            resTimeHist_Binned_PB = ResTimeHist_PB;
        end

        % Remove values < 0 from the histogram for fitting only
        idx = find(resTimeHist_Binned(:,2) < 0);
        resTimeHist_Binned(idx,2) = 0;
        idx = find(resTimeHist_Binned(:,2) <= 0);
        resTimeHist_Binned2 = resTimeHist_Binned;
        resTimeHist_Binned2(idx,:) = [];
        % Fit the binned residence time histogram with exponentials;
        [fitpar,espSigma, ~]= ExpDecay_fit_resTime(resTimeHist_Binned2, 1 );
        fit = [resTimeHist_Binned(:,1),ExpDecay_fun_resTime(fitpar,resTimeHist_Binned(:,1))];
        [fitpar2,espSigma2, ~]= ExpDecay_2Cmp_fit_resTime(resTimeHist_Binned2, [0.1 0.01]);

        Esp_Coef_1 = [fitpar2(1) fitpar2(4)*fitpar2(3)];
        Esp_Fit(:,1) = ExpDecay_fun_resTime(Esp_Coef_1,resTimeHist_Binned(:,1));

        Esp_Coef_2 = [fitpar2(2) fitpar2(4)*(1 - fitpar2(3))];
        Esp_Fit(:,2) = ExpDecay_fun_resTime(Esp_Coef_2,resTimeHist_Binned(:,1));

        fit2 = [resTimeHist_Binned(:,1),ExpDecay_2Cmp_fun_resTime(fitpar2,resTimeHist_Binned(:,1)),Esp_Fit];

        handles.Hist.Res(:, 1:2) = resTimeHist_Binned;
        handles.Hist.Res(:,3) = fit(:,2);
        handles.Hist.Res(:,4) = fit2(:,2);
        handles.Hist.Res(:,5) = fit2(:,3);
        handles.Hist.Res(:,6) = fit2(:,4);

        [~, pvalRes,FstatRes] = FtestModelCompare(resTimeHist_Binned(:,2),fit(:,2),fit2(:,2),2,4);

        handles.FitPar.Res = [fitpar, fitpar2, PartialBoundFraction, pvalRes; espSigma, espSigma2, BFerror,FstatRes];
        idx = find(resTimeHist_Binned_PB(:,2) < 0);
        resTimeHist_Binned_PB(idx,2) = 0;
        idx = find(resTimeHist_Binned_PB(:,2) <= 0);
        resTimeHist_Binned_PB2 = resTimeHist_Binned_PB;
        resTimeHist_Binned_PB2(idx,:) = [];
        % Fit the binned residence time histogram with exponentials;
        [fitpar_PB,espSigma_PB, fit]= ExpDecay_fit_resTime(resTimeHist_Binned_PB2, 1 );
        fit_PB = [resTimeHist_Binned(:,1),ExpDecay_fun_resTime(fitpar_PB,resTimeHist_Binned_PB(:,1))];

        [fitpar2_PB,espSigma2_PB, fit2]= ExpDecay_2Cmp_fit_resTime(resTimeHist_Binned_PB2, [0.1 0.01]);

        Esp_Coef_1 = [fitpar2_PB(1) fitpar2_PB(4)*fitpar2_PB(3)];
        Esp_Fit(:,1) = ExpDecay_fun_resTime(Esp_Coef_1,resTimeHist_Binned_PB(:,1));

        Esp_Coef_2 = [fitpar2_PB(2) fitpar2_PB(4)*(1 - fitpar2_PB(3))];
        Esp_Fit(:,2) = ExpDecay_fun_resTime(Esp_Coef_2,resTimeHist_Binned_PB(:,1));

        fit2_PB = [resTimeHist_Binned_PB(:,1),ExpDecay_2Cmp_fun_resTime(fitpar2_PB,resTimeHist_Binned_PB(:,1)),Esp_Fit];
        handles.Hist.Res_PB(:, 1:2) = resTimeHist_Binned_PB;
        handles.Hist.Res_PB(:,3) = fit_PB(:,2);
        handles.Hist.Res_PB(:,4) = fit2_PB(:,2);
        handles.Hist.Res_PB(:,5) = fit2_PB(:,3);
        handles.Hist.Res_PB(:,6) = fit2_PB(:,4);
        [~, pvalRes_PB,FstatRes_PB] = FtestModelCompare(resTimeHist_Binned_PB(:,2),fit_PB(:,2),fit2_PB(:,2),2,4);
        handles.FitPar.Res_PB = [fitpar_PB, fitpar2_PB, PartialBoundFraction, pvalRes_PB; espSigma_PB, espSigma2_PB, BFerror, FstatRes_PB];


        HistAxes = findobj ('Tag', 'hist_axes');
        axes(HistAxes);
        if strcmp(answer,'Yes');
            resTimeHist_Binned = resTimeHist_Binned_PB;
            fit = fit_PB;
            fit2 = fit2_PB;
            fitpar = fitpar_PB;
            fitpar2 = fitpar2_PB;
            espSigma = espSigma_PB;
            espSigma2 = espSigma2_PB;
            pvalRes = pvalRes_PB;
            FstatRes = FstatRes_PB;
        end
        if handles.Flag == 1
            bar(resTimeHist_Binned(:,1), resTimeHist_Binned(:,2),'EdgeColor', [0 0 0], ...
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
                ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction,3), ' \pm ',...
                num2str(BFerror,2)],...
                ['One Component fit k_{off} = ',num2str(fitpar(1),3),...
                ' \pm ', num2str(espSigma(1),2), ' s^{-1}; C_{eq} = ', num2str(fitpar(2),2)],...
                ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
                ' \pm ', num2str(espSigma2(1),2), ' s^{-1}'], ...
                ['k_{off, 2} = ', num2str(fitpar2(2),3),...
                ' \pm ', num2str(espSigma2(2),2), ' s^{-1}; f_1 = ', num2str(fitpar2(3),3),...
                '\pm', num2str(espSigma2(3),2), ' C_{eq} = ',  num2str(fitpar2(4),2)]});
            legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit','Component 1', 'Component 2');
            box on;
            set(handles.FstatText,'String',['F stat: ', num2str(FstatRes,4)]);
            set(handles.pValText,'String',['p-value: ', num2str(pvalRes,4)]);

        end




    % end

    % if handles.Flag == 2;
    %     set(handles.What2Plot,'Value',2);
        % Logarithmic sampling of the survival prob
        if handles.Parameters.bin ~= 0

    %         if strcmp(answer,'Yes')
                TimePoints_PB = logspace(log10(min(CumTrackLengthHist_PB(:,1))),...
                    log10(max(CumTrackLengthHist_PB(:,1))),handles.Parameters.bin);
    %         else
                TimePoints = logspace(log10(min(CumTrackLengthHist(:,1))),...
                    log10(max(CumTrackLengthHist(:,1))),handles.Parameters.bin);
    %         end
            TrLog =[];

            TimePoints = round(TimePoints*1000)/1000;
            TimePoints_PB = round(TimePoints_PB*1000)/1000;
            for i = 1:length(TimePoints)

                idx = find(CumTrackLengthHist(:,1) <= TimePoints(i)+0.0001,1,'last');
                TrLog(i,:) = CumTrackLengthHist(idx,:);
            end
            for i = 1:length(TimePoints_PB)

                idx = find(CumTrackLengthHist_PB(:,1) <= TimePoints_PB(i)+0.0001,1,'last');
                TrLog_PB(i,:) = CumTrackLengthHist_PB(idx,:);
            end

        else
            TrLog = CumTrackLengthHist;
            TrLog_PB = CumTrackLengthHist_PB;
        end

        tpoint1 = TrLog(1,1);
        TrLog(:,1) = TrLog(:,1) - tpoint1;


    %     ind2fit = 1:length(TrLog);
        % Fit exponentials to the Survival probability
        ind2fit = 1:size(TrLog,1);
        if handles.Parameters.bin == 0

            [fitpar,espSigma, ~]= ExpDecay_fit(TrLog(ind2fit,:), 1/(tpoint1*10));
            fit = [TrLog(:,1),ExpDecay_fun(fitpar,TrLog(:,1))];
            [fitpar2,espSigma2, ~]= ExpDecay_2Cmp_fit(TrLog(ind2fit,:), [1/(tpoint1*10) 1/(tpoint1*100)]);
            Esp_Coef_1 = [fitpar2(1) fitpar2(4)*fitpar2(3)];
            Esp_Fit(:,1) = ExpDecay_fun(Esp_Coef_1,TrLog(:,1));

            Esp_Coef_2 = [fitpar2(2) fitpar2(4)*(1 - fitpar2(3))];
            Esp_Fit(:,2) = ExpDecay_fun(Esp_Coef_2,TrLog(:,1));
            fit2 = [TrLog(:,1),ExpDecay_2Cmp_fun(fitpar2,TrLog(:,1)),Esp_Fit];
        else
            [fitpar,espSigma, fit]= ExpDecay_fit(TrLog(ind2fit,:), 1/(tpoint1*10));
            [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit(TrLog(ind2fit,:), [1/(tpoint1*10) 1/(tpoint1*100)]);
        end

        fit(:,1) = fit(:,1) + tpoint1;
        fit2(:,1) = fit2(:,1) + tpoint1;
        TrLog(:,1) = TrLog(:,1) + tpoint1;

        tpoint1 = TrLog_PB(1,1);
        TrLog_PB(:,1) = TrLog_PB(:,1) - tpoint1;



        % Fit exponentials to the Survival probability
        if handles.Parameters.bin == 0
            [fitpar_PB,espSigma_PB, ~]= ExpDecay_fit(TrLog_PB(ind2fit,:), 1/(tpoint1*10));
            fit_PB = [TrLog_PB(:,1),ExpDecay_fun(fitpar_PB,TrLog_PB(:,1))];
            [fitpar2_PB,espSigma2_PB, ~]= ExpDecay_2Cmp_fit(TrLog_PB(ind2fit,:), [1/(tpoint1*10) 1/(tpoint1*100)]);
            Esp_Coef_1 = [fitpar2_PB(1) fitpar2_PB(4)*fitpar2_PB(3)];
            Esp_Fit(:,1) = ExpDecay_fun(Esp_Coef_1,TrLog_PB(:,1));

            Esp_Coef_2 = [fitpar2_PB(2) fitpar2_PB(4)*(1 - fitpar2_PB(3))];
            Esp_Fit(:,2) = ExpDecay_fun(Esp_Coef_2,TrLog_PB(:,1));
            fit2_PB = [TrLog_PB(:,1),ExpDecay_2Cmp_fun(fitpar2_PB,TrLog_PB(:,1)),Esp_Fit];
        else
            [fitpar_PB,espSigma_PB, fit_PB]= ExpDecay_fit(TrLog_PB(ind2fit,:), 1/(tpoint1*10));
            [fitpar2_PB,espSigma2_PB, fit2_PB]= ExpDecay_2Cmp_fit(TrLog_PB(ind2fit,:), [1/(tpoint1*10) 1/(tpoint1*100)]);
        end
        fit_PB(:,1) = fit_PB(:,1) + tpoint1;
        fit2_PB(:,1) = fit2_PB(:,1) + tpoint1;
        TrLog_PB(:,1) = TrLog_PB(:,1) + tpoint1;

        handles.Hist.Surv(:, 1:2) = TrLog;
        handles.Hist.Surv(:,3) = fit(:,2);
        handles.Hist.Surv(:,4) = fit2(:,2);
        handles.Hist.Surv(:,5) = fit2(:,3);
        handles.Hist.Surv(:,6) = fit2(:,4);
        [~, pvalSurv, FstatSurv] = FtestModelCompare(TrLog(:,2),fit(:,2),fit2(:,2),2,4);
        handles.FitPar.Surv = [fitpar, fitpar2, PartialBoundFraction, pvalSurv; espSigma, espSigma2, BFerror, FstatSurv];

        handles.Hist.Surv_PB(:, 1:2) = TrLog_PB;
        handles.Hist.Surv_PB(:,3) = fit_PB(:,2);
        handles.Hist.Surv_PB(:,4) = fit2_PB(:,2);
        handles.Hist.Surv_PB(:,5) = fit2_PB(:,3);
        handles.Hist.Surv_PB(:,6) = fit2_PB(:,4);
        [~, pvalSurv_PB, FstatSurv_PB] = FtestModelCompare(TrLog_PB(:,2),fit_PB(:,2),fit2_PB(:,2),2,4);
        handles.FitPar.Surv_PB = [fitpar_PB, fitpar2_PB, PartialBoundFraction, pvalSurv_PB; espSigma_PB, espSigma2_PB, BFerror, FstatSurv_PB];
        if strcmp(answer,'Yes');
            TrLog = TrLog_PB;
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
            HistAxes = findobj ('Tag', 'hist_axes');
            axes(HistAxes);
            plot(TrLog(:,1), TrLog(:,2),'ok');
            hold on;
            plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
            plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
            plot(fit2(:,1),fit2(:,3),'b', 'LineWidth', 1);
            plot(fit2(:,1),fit2(:,4),'b--', 'LineWidth', 1);
            hold off;
            xlabel('Time [s]');
            ylabel('Counts');
            title({['Fit of the survival probability',...
                ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction,3), ' \pm ',...
                num2str(BFerror,2)],...
                ['One Component fit k_{off} = ',num2str(fitpar(1),3),...
                ' \pm ', num2str(espSigma(1),2), ' s^{-1}; C_{eq} = ', num2str(fitpar(2),2)],...
                ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
                ' \pm ', num2str(espSigma2(1),2), ' s^{-1}'], ...
                ['k_{off, 2} = ', num2str(fitpar2(2),3),...
                ' \pm ', num2str(espSigma2(2),2), ' s^{-1}; f_1 = ', num2str(fitpar2(3),3),...
                '\pm', num2str(espSigma2(3),2), 'C_{eq} = ',  num2str(fitpar2(4),2)]});
            legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit');
            box on;
            set(handles.FstatText,'String',['F stat: ', num2str(FstatSurv,4)]);
            set(handles.pValText,'String',['p-value: ', num2str(pvalSurv,4)]);
        end





    set(handles.totalTracksText,'String',['Total Tracks: ',num2str(handles.totalTracks)]);
    set(handles.boundParticlesText,'String',['Bound Particles: ',num2str(handles.BoundMolecules)]);
    set(handles.totalParticlesText,'String',['Total Particles: ',num2str(handles.TotalMolecules)]);
    set(handles.totalROIsText,'String',['Total ROIs: ', num2str(handles.nROIs)]);
    

    figname = get(gcf,'Name');
    figname = [figname ' ' ROIChoice{1}];
    set(gcf,'Name',figname);
    set(handles.PieChart,'Enable','on');
  
    set(handles.frame2frameMaxText,'String',['Max. Frame to Frame Jump: ', num2str(handles.Parameters.ThreshL), ' um']);
    set(handles.End2EndMaxText,'String',['Max. End to End Jump: ', num2str(handles.Parameters.ThreshH), ' um']);
    set(handles.NminText,'String',['Min. Bound Frames (Nmin): ', num2str(handles.Parameters.minBoundFrames)]);
end

    % Update handles structure
    guidata(hObject, handles);
t_elap = toc;
disp(['Execution Time = ' num2str(t_elap), ' s']);

% UIWAIT makes AnalyzeResTimes_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AnalyzeResTimes_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in CopyHist.
function CopyHist_Callback(hObject, eventdata, handles)
% hObject    handle to CopyHist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.What2Plot,'Value') == 1 && ~get(handles.PBCorrect,'Value')
    Hist = handles.Hist.Res;
elseif get(handles.What2Plot,'Value') == 1 && get(handles.PBCorrect,'Value')
    Hist = handles.Hist.Res_PB;
elseif get(handles.What2Plot,'Value') == 2 && ~get(handles.PBCorrect,'Value')
    Hist = handles.Hist.Surv;
elseif get(handles.What2Plot,'Value') == 2 && get(handles.PBCorrect,'Value')
    Hist = handles.Hist.Surv_PB;
end
num2clip(Hist);


% --- Executes on button press in CopyAnPar.
function CopyAnPar_Callback(hObject, eventdata, handles)
% hObject    handle to CopyAnPar (see GCBO)
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
    FitPar = handles.FitPar.Res;
elseif get(handles.What2Plot,'Value') == 1 && get(handles.PBCorrect,'Value')
    FitPar = handles.FitPar.Res_PB;
elseif get(handles.What2Plot,'Value') == 2 && ~get(handles.PBCorrect,'Value')
    FitPar = handles.FitPar.Surv;
elseif get(handles.What2Plot,'Value') == 2 && get(handles.PBCorrect,'Value')
    FitPar = handles.FitPar.Surv_PB;
end

num2clip(FitPar);




% --- Executes on button press in SaveFigure.
function SaveFigure_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Prepare data for figure;

CumNParticles = handles.NParticles;
BleachRates = handles.BleachRates;
CumTrackLengthHist = handles.CumHist;
if get(handles.What2Plot,'Value') == 1 && ~get(handles.PBCorrect,'Value')
    Hist = handles.Hist.Res;
    FitParAll = handles.FitPar.Res;
elseif get(handles.What2Plot,'Value') == 1 && get(handles.PBCorrect,'Value')
    Hist = handles.Hist.Res_PB;
    FitParAll = handles.FitPar.Res_PB;
elseif get(handles.What2Plot,'Value') == 2 && ~get(handles.PBCorrect,'Value')
    Hist = handles.Hist.Surv;
    FitParAll = handles.FitPar.Surv;
elseif get(handles.What2Plot,'Value') == 2 && get(handles.PBCorrect,'Value')
    Hist = handles.Hist.Surv_PB;
    FitParAll = handles.FitPar.Surv_PB;
end


% Hist = handles.Hist;
Flag = get(handles.What2Plot,'Value');

fitpar = FitParAll(1, 1:2);
fitpar2 = FitParAll(1, 3:6);
PartialBoundFraction = FitParAll(1,7);

espSigma = FitParAll(2, 1:2);
espSigma2 = FitParAll(2, 3:6);
BFerror = FitParAll(2, 7);




% Plot the photobleaching curve
scrsz = get(0,'ScreenSize');
h1 = figure('Position',[.15*scrsz(3) .15*scrsz(4) .5*scrsz(3) .4*scrsz(4)]);



subplot(1,2,1);
plot(CumNParticles(:,1), CumNParticles(:,2)/max(CumNParticles(:,4)), '+', 'MarkerEdgeColor', [0.5 0.5 0.5]');
hold on
plot(CumTrackLengthHist(:,1), CumTrackLengthHist(:,2)/max(CumTrackLengthHist(:,2)),'ok');
hold off;
xlabel('Time [s]');
ylabel('Normalized Counts');
legend('Photobleaching Decay', 'Bound molecules Decay');
hold on;
plot(CumNParticles(:,3), CumNParticles(:,4)/max(CumNParticles(:,4)),'r');

title({'Photobleaching kinetics:', ['k_{b1} = ', num2str(BleachRates(1), 3),...
    ' s^{-1}, k_{b2} = ', num2str(BleachRates(2), 3), ' s^{-1}, f_1 = ',  ...
    num2str(BleachRates(3), 3)]})

hold off;
box on;

% Plot the residence time/ cumulative histogram.

subplot(1,2,2);
plot(Hist(:,1), Hist(:,2),'ok');
hold on;
plot(Hist(:,1),Hist(:,3),'g', 'LineWidth', 1);
plot(Hist(:,1),Hist(:,4),'r', 'LineWidth', 1);
plot(Hist(:,1),Hist(:,5),'b', 'LineWidth', 1);
plot(Hist(:,1),Hist(:,6),'b--', 'LineWidth', 1);
hold off;
box on;
if Flag == 1
    xlabel('Residence Time [s]')
    ylabel('Counts');
    title({['Fit of the residence time histogram',...
        ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction,3), ' \pm ',...
        num2str(BFerror,2)],...
        ['One Component fit k_{off} = ',num2str(fitpar(1),3),...
        ' \pm ', num2str(espSigma(1),2), ' s; C_eq = ', num2str(fitpar(2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
        ' \pm ', num2str(espSigma2(1),2), 's'], ...
        ['k_{off, 2} = ', num2str(fitpar2(2),3),...
        ' \pm ', num2str(espSigma2(2),2), 's; f_1 = ', num2str(fitpar2(3),3),...
        '\pm', num2str(espSigma2(3),2), 'C_eq = ',  num2str(fitpar2(4),2)]});
    legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit','Component 1','Component 2');
    box on;
    DefaultName = handles.PathName;
    DefaultName = [DefaultName, 'GlobalResTHist'];
    
elseif Flag == 2;
    xlabel('Time [s]');
    ylabel('Counts');
    title({['Fit of the survival probability',...
        ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction,3), ' \pm ',...
        num2str(BFerror,2)],...
        ['One Component fit k_{off} = ',num2str(fitpar(1),3),...
        ' \pm ', num2str(espSigma(1),2), ' s; C_eq = ', num2str(fitpar(2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
        ' \pm ', num2str(espSigma2(1),2), 's'], ...
        ['k_{off, 2} = ', num2str(fitpar2(2),3),...
        ' \pm ', num2str(espSigma2(2),2), 's; f_1 = ', num2str(fitpar2(3),3),...
        '\pm', num2str(espSigma2(3),2), 'C_eq = ',  num2str(fitpar2(4),2)]});
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


% --- Executes on selection change in What2Plot.
function What2Plot_Callback(hObject, eventdata, handles)
% hObject    handle to What2Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns What2Plot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from What2Plot
PartialBoundFraction = handles.PartBoundFrac;
BFerror = handles.BFerror;
if get(hObject,'Value') == 1;
    

    if get(handles.PBCorrect,'Value') == 1 
        Hist = handles.Hist.Res_PB;
        pval = handles.FitPar.Res_PB(1,8);
        Fstat = handles.FitPar.Res_PB(2,8);
        
    else
        Hist = handles.Hist.Res;
        pval = handles.FitPar.Res(1,8);
        Fstat = handles.FitPar.Res(2,8);
    end
    
    if size(Hist,2) < 3

        % Bin the histogram according to the selected settings;

        binFactor = handles.Parameters.bin + 1;
        nBins = floor(length(Hist(:,1))/binFactor);
        
            TimeBins = reshape(Hist(1:nBins*binFactor,1), binFactor, nBins);
            CountsBins = reshape(Hist(1:nBins*binFactor,2), binFactor, nBins);
        

        if binFactor ~= 1
            resTimeHist_Binned = [];
            resTimeHist_Binned(:,1) = mean(TimeBins);
            resTimeHist_Binned(:,2) = sum(CountsBins);
        else
            resTimeHist_Binned = Hist;
        end

        % Remove values < 0 from the histogram
        idx = find(resTimeHist_Binned(:,2) <= 0);
        resTimeHist_Binned(idx,:) = [];
        % Fit the binned residence time histogram with exponentials;
        [fitpar,espSigma, fit]= ExpDecay_fit_resTime(resTimeHist_Binned, 1 );
        [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit_resTime(resTimeHist_Binned, [1 0.1]);
        if get(handles.PBCorrect,'Value') == 1 
            handles.Hist.Res_PB(:, 1:2) = resTimeHist_Binned;
            handles.Hist.Res_PB(:,3) = fit(:,2);
            handles.Hist.Res_PB(:,4) = fit2(:,2);
            handles.Hist.Res_PB(:,5) = fit2(:,3);
            handles.Hist.Res_PB(:,6) = fit2(:,4);
            
            handles.FitPar.Res_PB = [fitpar, fitpar2, PartialBoundFraction; espSigma, espSigma2, BFerror];
        else
            handles.Hist.Res(:, 1:2) = resTimeHist_Binned;
            handles.Hist.Res(:,3) = fit(:,2);
            handles.Hist.Res(:,4) = fit2(:,2);
            handles.Hist.Res(:,5) = fit2(:,3);
            handles.Hist.Res(:,6) = fit2(:,4);
            
            handles.FitPar.Res = [fitpar, fitpar2, PartialBoundFraction; espSigma, espSigma2, BFerror];
        end
        
    else
        resTimeHist_Binned = Hist;
        fit = [Hist(:,1), Hist(:,3)];
        fit2 = [Hist(:,1),Hist(:,4:6)];
        if get(handles.PBCorrect,'Value') == 1 
            fitpar = handles.FitPar.Res_PB(1,1:2);
            espSigma = handles.FitPar.Res_PB(2,1:2);
            fitpar2 = handles.FitPar.Res_PB(1,3:6);
            espSigma2 = handles.FitPar.Res_PB(2,3:6);
        else
            fitpar = handles.FitPar.Res(1,1:2);
            espSigma = handles.FitPar.Res(2,1:2);
            fitpar2 = handles.FitPar.Res(1,3:6);
            espSigma2 = handles.FitPar.Res(2,3:6);
        end
    end
        
    
    
    
%     HistAxes = findobj ('Tag', 'hist_axes');
    axes(handles.hist_axes);
    
    bar(resTimeHist_Binned(:,1), resTimeHist_Binned(:,2),'EdgeColor', [0 0 0], ...
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
        ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction,3), ' \pm ',...
        num2str(BFerror,2)],...
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
        Hist = handles.Hist.Surv_PB;
        pval = handles.FitPar.Surv_PB(1,8);
        Fstat = handles.FitPar.Surv_PB(2,8);
    else
        Hist = handles.Hist.Surv;
        pval = handles.FitPar.Surv(1,8);
        Fstat = handles.FitPar.Surv(2,8);
    end
    
    if size(Hist,2) < 3
    
        % Logarithmic sampling of the survival prob
        if handles.Parameters.bin ~= 0
            
            
            TimePoints = logspace(log10(min(Hist(:,1))),...
                log10(max(Hist(:,1))),handles.Parameters.bin);
            
            TrLog =[];
            
            TimePoints = round(TimePoints*1000)/1000;
            for i = 1:length(TimePoints)
                
                idx = find(Hist(:,1) <= TimePoints(i)+0.0001,1,'last');
                TrLog(i,:) = Hist(idx,:);
            end
            
        else
            TrLog = Hist;
        end
        
        tpoint1 = TrLog(1,1);
        TrLog(:,1) = TrLog(:,1) - tpoint1;
        
        
        
        
        % Fit exponentials to the Survival probability
        [fitpar,espSigma, fit]= ExpDecay_fit(TrLog, 1/(tpoint1*10));
        [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit(TrLog, [1/(tpoint1*10) 1/(tpoint1*100)]);
        fit(:,1) = fit(:,1) + tpoint1;
        fit2(:,1) = fit2(:,1) + tpoint1;
        TrLog(:,1) = TrLog(:,1) + tpoint1;
        if get(handles.PBCorrect,'Value') == 0
            handles.Hist.Surv(:, 1:2) = TrLog;
            handles.Hist.Surv(:,3) = fit(:,2);
            handles.Hist.Surv(:,4) = fit2(:,2);
            handles.Hist.Surv(:,5) = fit2(:,3);
            handles.Hist.Surv(:,6) = fit2(:,4);
            handles.FitPar.Surv = [fitpar, fitpar2, PartialBoundFraction; espSigma, espSigma2, BFerror];
        else
            handles.Hist.Surv_PB(:, 1:2) = TrLog;
            handles.Hist.Surv_PB(:,3) = fit(:,2);
            handles.Hist.Surv_PB(:,4) = fit2(:,2);
            handles.Hist.Surv_PB(:,5) = fit2(:,3);
            handles.Hist.Surv_PB(:,6) = fit2(:,4);
            handles.FitPar.Surv_PB = [fitpar, fitpar2, PartialBoundFraction; espSigma, espSigma2, BFerror];
        end
    else
        TrLog = Hist;
        fit = [Hist(:,1),Hist(:,3)];
        fit2 = [Hist(:,1),Hist(:,4:6)];
        if get(handles.PBCorrect,'Value') == 1
            fitpar = handles.FitPar.Surv_PB(1,1:2);
            espSigma = handles.FitPar.Surv_PB(2,1:2);
            fitpar2 = handles.FitPar.Surv_PB(1,3:6);
            espSigma2 = handles.FitPar.Surv_PB(2,3:6);
        else
            fitpar = handles.FitPar.Surv(1,1:2);
            espSigma = handles.FitPar.Surv(2,1:2);
            fitpar2 = handles.FitPar.Surv(1,3:6);
            espSigma2 = handles.FitPar.Surv(2,3:6);
        end
    end
    
%     HistAxes = findobj ('Tag', '');
    axes(handles.hist_axes);
    plot(TrLog(:,1), TrLog(:,2),'ok');
    hold on;
    plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,3),'b', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,4),'b--', 'LineWidth', 1);
    hold off;
    xlabel('Time [s]');
    ylabel('Counts');
    title({['Fit of the survival probability',...
        ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction,3), ' \pm ',...
        num2str(BFerror,2)],...
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
handles.CumHist = Hist;
set(handles.FstatText,'String',['F stat: ', num2str(Fstat,4)]);
set(handles.pValText,'String',['p-value: ', num2str(pval,4)]);
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
PartialBoundFraction = handles.PartBoundFrac;
BFerror = handles.BFerror;
if get(handles.What2Plot,'Value') == 1;
    

    if get(handles.PBCorrect,'Value') == 1 
        Hist = handles.Hist.Res_PB;
        pval = handles.FitPar.Res(1,8);
        Fstat = handles.FitPar.Res(2,8);
    else
        Hist = handles.Hist.Res;
        pval = handles.FitPar.Res_PB(1,8);
        Fstat = handles.FitPar.Res_PB(2,8);
    end
    
    if size(Hist,2) < 3

        % Bin the histogram according to the selected settings;

        binFactor = handles.Parameters.bin + 1;
        nBins = floor(length(ResTimeHist(:,1))/binFactor);
        
            TimeBins = reshape(Hist(1:nBins*binFactor,1), binFactor, nBins);
            CountsBins = reshape(Hist(1:nBins*binFactor,2), binFactor, nBins);
        

        if binFactor ~= 1
            resTimeHist_Binned = [];
            resTimeHist_Binned(:,1) = mean(TimeBins);
            resTimeHist_Binned(:,2) = sum(CountsBins);
        else
            resTimeHist_Binned = ResTimeHist;
        end

        % Remove values < 0 from the histogram
        idx = find(resTimeHist_Binned(:,2) <= 0);
        resTimeHist_Binned(idx,:) = [];
        % Fit the binned residence time histogram with exponentials;
        [fitpar,espSigma, fit]= ExpDecay_fit_resTime(resTimeHist_Binned, 1 );
        [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit_resTime(resTimeHist_Binned, [1 0.1]);
        if get(handles.PBCorrect,'Value') == 1 
            handles.Hist.Res_PB(:, 1:2) = resTimeHist_Binned;
            handles.Hist.Res_PB(:,3) = fit(:,2);
            handles.Hist.Res_PB(:,4) = fit2(:,2);
            handles.Hist.Res_PB(:,5) = fit2(:,3);
            handles.Hist.Res_PB(:,6) = fit2(:,4);
            
            handles.FitPar.Res_PB = [fitpar, fitpar2, PartialBoundFraction; espSigma, espSigma2, BFerror];
        else
            handles.Hist.Res(:, 1:2) = resTimeHist_Binned;
            handles.Hist.Res(:,3) = fit(:,2);
            handles.Hist.Res(:,4) = fit2(:,2);
            handles.Hist.Res(:,5) = fit2(:,3);
            handles.Hist.Res(:,6) = fit2(:,4);
            
            handles.FitPar.Res = [fitpar, fitpar2, PartialBoundFraction; espSigma, espSigma2, BFerror];
        end
        
    else
        resTimeHist_Binned = Hist;
        fit = [Hist(:,1),Hist(:,3)];
        fit2 = [Hist(:,1),Hist(:,4:6)];
        if get(handles.PBCorrect,'Value') == 1 
            fitpar = handles.FitPar.Res_PB(1,1:2);
            espSigma = handles.FitPar.Res_PB(2,1:2);
            fitpar2 = handles.FitPar.Res_PB(1,3:6);
            espSigma2 = handles.FitPar.Res_PB(2,3:6);
        else
            fitpar = handles.FitPar.Res(1,1:2);
            espSigma = handles.FitPar.Res(2,1:2);
            fitpar2 = handles.FitPar.Res(1,3:6);
            espSigma2 = handles.FitPar.Res(2,3:6);
        end
    end
        
    
    
    
%     HistAxes = findobj ('Tag', 'hist_axes');
    axes(handles.hist_axes);
    bar(resTimeHist_Binned(:,1), resTimeHist_Binned(:,2),'EdgeColor', [0 0 0], ...
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
        ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction,3), ' \pm ',...
        num2str(BFerror,2)],...
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
        Hist = handles.Hist.Surv_PB;
        pval = handles.FitPar.Surv_PB(1,8);
        Fstat = handles.FitPar.Surv_PB(2,8);
    else
        Hist = handles.Hist.Surv;
        pval = handles.FitPar.Surv(1,8);
        Fstat = handles.FitPar.Surv(2,8);
    end
    
    if size(Hist,2) < 3
    
        % Logarithmic sampling of the survival prob
        if handles.Parameters.bin ~= 0
            
            
            TimePoints = logspace(log10(min(Hist(:,1))),...
                log10(max(Hist(:,1))),handles.Parameters.bin);
            
            TrLog =[];
            
            TimePoints = round(TimePoints*1000)/1000;
            for i = 1:length(TimePoints)
                
                idx = find(Hist(:,1) <= TimePoints(i)+0.0001,1,'last');
                TrLog(i,:) = Hist(idx,:);
            end
            
        else
            TrLog = Hist;
        end
        
        tpoint1 = TrLog(1,1);
        TrLog(:,1) = TrLog(:,1) - tpoint1;
        
        
        
        
        % Fit exponentials to the Survival probability
        [fitpar,espSigma, fit]= ExpDecay_fit(TrLog, 1/(tpoint1*10));
        [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit(TrLog, [1/(tpoint1*10) 1/(tpoint1*100)]);
        fit(:,1) = fit(:,1) + tpoint1;
        fit2(:,1) = fit2(:,1) + tpoint1;
        TrLog(:,1) = TrLog(:,1) + tpoint1;
        if get(handles.PBCorrect,'Value') == 0
            handles.Hist.Surv(:, 1:2) = TrLog;
            handles.Hist.Surv(:,3) = fit(:,2);
            handles.Hist.Surv(:,4) = fit2(:,2);
            handles.Hist.Surv(:,5) = fit2(:,3);
            handles.Hist.Surv(:,6) = fit2(:,4);
            handles.FitPar.Surv = [fitpar, fitpar2, PartialBoundFraction; espSigma, espSigma2, BFerror];
        else
            handles.Hist.Surv_PB(:, 1:2) = TrLog;
            handles.Hist.Surv_PB(:,3) = fit(:,2);
            handles.Hist.Surv_PB(:,4) = fit2(:,2);
            handles.Hist.Surv_PB(:,5) = fit2(:,3);
            handles.Hist.Surv_PB(:,6) = fit2(:,4);
            handles.FitPar.Surv_PB = [fitpar, fitpar2, PartialBoundFraction; espSigma, espSigma2, BFerror];
        end
    else
        TrLog = Hist;
        fit = [Hist(:,1),Hist(:,3)];
        fit2 = [Hist(:,1),Hist(:,4:6)];
        if get(handles.PBCorrect,'Value') == 1
            fitpar = handles.FitPar.Surv_PB(1,1:2);
            espSigma = handles.FitPar.Surv_PB(2,1:2);
            fitpar2 = handles.FitPar.Surv_PB(1,3:6);
            espSigma2 = handles.FitPar.Surv_PB(2,3:6);
        else
            fitpar = handles.FitPar.Surv(1,1:2);
            espSigma = handles.FitPar.Surv(2,1:2);
            fitpar2 = handles.FitPar.Surv(1,3:6);
            espSigma2 = handles.FitPar.Surv(2,3:6);
        end
    end
    
%     HistAxes = findobj ('Tag', 'hist_axes');
    axes(handles.hist_axes);
    plot(TrLog(:,1), TrLog(:,2),'ok');
    hold on;
    plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,3),'b', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,4),'b--', 'LineWidth', 1);
    hold off
    xlabel('Time [s]');
    ylabel('Counts');
    title({['Fit of the survival probability',...
        ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction,3), ' \pm ',...
        num2str(BFerror,2)],...
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
handles.CumHist = Hist;
set(handles.FstatText,'String',['F stat: ', num2str(Fstat,4)]);
set(handles.pValText,'String',['p-value: ', num2str(pval,4)]);
guidata(hObject,handles);


% --- Executes on button press in BoxPlotButton.
function BoxPlotButton_Callback(hObject, eventdata, handles)
% hObject    handle to BoxPlotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Data = handles.Hist.Res_PB(:,1:2);
SurvData = handles.Hist.Surv_PB(:,1:2);

ResVec = Data(:,2);
ResVec = Data(:,2)./sum(Data(:,2));
% ResVec = ResVec.*SurvData(1,2);

% ResVecN = ResVec.*handles.BoundMolecules;
ResVecN = ResVec.*handles.totalTracks;

ResVecN = round(ResVecN);
ResTime = Data(:,1);


    

BoxPlotData = [];
if handles.FitPar.Surv_PB(1,8) > 0.01
    
    for i = 1:length(ResTime)
        vecTmp = [ResTime(i)*ones(ResVecN(i),1), ones(ResVecN(i),1)];
        BoxPlotData = [BoxPlotData; vecTmp];
        
    end
else
    
    f1 = handles.FitPar.Surv_PB(1,5);
    f2 = (1 - f1);
    a = handles.FitPar.Surv_PB(1,6);
    
    cutoffInd = find(SurvData(:,2) >= a*f2,1,'last');
    for i = 1:cutoffInd
        vecTmp = [ResTime(i)*ones(ResVecN(i),1), ones(ResVecN(i),1)];
        BoxPlotData = [BoxPlotData; vecTmp];
        
    end
    
    for i = cutoffInd+1:length(ResTime)
        vecTmp = [ResTime(i)*ones(ResVecN(i),1), 2*ones(ResVecN(i),1)];
        BoxPlotData = [BoxPlotData; vecTmp];
        
    end
end
figure; boxplot(BoxPlotData(:,1),BoxPlotData(:,2));
num2clip(BoxPlotData);
msgbox('Box Plot data saved to clipboard','Box plot data');


% --------------------------------------------------------------------
function Untitled1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SaveData_Callback(hObject, eventdata, handles)
% hObject    handle to SaveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

defPath = handles.PathName;

[fname,pname] = uiputfile('*.mat','Save Residence Time Analysis',[defPath filesep 'ResTimeAnalysis.mat']);

if fname ~= 0
    Results.Hist = handles.Hist;
    Results.FitPar = handles.FitPar;
    Results.NParticles = handles.NParticles;
    Results.BleachRates = handles.BleachRates;
    Results.Parameters = handles.Parameters;
    Results.totalTracks = handles.totalTracks;
    Results.BoundMolecules = handles.BoundMolecules;
    Results.TotalMolecules = handles.TotalMolecules;
    Results.totalROIs = handles.nROIs;
    Results.PartBoundFrac = handles.PartBoundFrac;
    save([pname,fname],'Results');
end


% --------------------------------------------------------------------
function LoadResData_Callback(hObject, eventdata, handles)
% hObject    handle to LoadResData (see GCBO)
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
    handles.NParticles = IN.Results.NParticles;
    handles.BleachRates = IN.Results.BleachRates;
    handles.Parameters = IN.Results.Parameters;
    handles.totalTracks = IN.Results.totalTracks;
    handles.BoundMolecules = IN.Results.BoundMolecules;
    handles.TotalMolecules = IN.Results.TotalMolecules;
    if isfield(IN.Results,'PartBoundFrac')
        handles.PartBoundFrac = IN.Results.PartBoundFrac;
    else
        handles.PartBoundFrac = IN.Results.BoundMolecules./IN.Results.TotalMolecules;
    end
    handles.CumHist = handles.Hist.Surv_PB;
    if isfield(IN.Results,'totalROIs');
        handles.nROIs = IN.Results.totalROIs;
        set(handles.totalROIsText,'String',['Total ROIs: ',num2str(handles.nROIs)]);
    else
        set(handles.totalROIsText,'String','Total ROIs: Not Available');
    end
    %Display Survival plot with PB correction by default
    set(handles.What2Plot,'Value',2);
    set(handles.PBCorrect,'Value',1);
    %Display track & particle counts
    set(handles.totalTracksText,'String',['Total Tracks: ', num2str(handles.totalTracks)]);
    set(handles.boundParticlesText,'String',['Bound Particles: ', num2str(handles.BoundMolecules)]);
    set(handles.totalParticlesText,'String',['Total Particles: ', num2str(handles.TotalMolecules)]);
    %Display F stat info
    set(handles.FstatText,'String',['F stat: ', num2str(handles.FitPar.Surv_PB(2,8))]);
    set(handles.pValText,'String', ['p-value: ', num2str(handles.FitPar.Surv_PB(1,8))]);
    
    %Plot bleaching curve
    axes(handles.photob_axes)
    plot(handles.NParticles(:,1), handles.NParticles(:,2)/max(handles.NParticles(:,4)), '+', 'MarkerEdgeColor', [0.5 0.5 0.5]');
    hold on
    plot(handles.Hist.Surv(:,1), handles.Hist.Surv(:,2)/max(handles.Hist.Surv(:,2)),'ok');
    plot(handles.NParticles(:,3), handles.NParticles(:,4)/max(handles.NParticles(:,4)),'r');
    hold off;
    
    box on;
    
    title({'Photobleaching kinetics:', ['k_{b1} = ', num2str(handles.BleachRates(1), 3),...
        ' s^{-1}, k_{b2} = ', num2str(handles.BleachRates(2), 3), ' s^{-1}, f_1 = ',  ...
        num2str(handles.BleachRates(3), 3)]})
    
    hold off;
    legend('Photobleaching Decay', 'Bound molecules Decay');
    
    xlabel('Time [s]');
    ylabel('Normalized Counts');
    
    %Plot the survival distribution and fits
    axes(handles.hist_axes)
    plot(handles.Hist.Surv_PB(:,1), handles.Hist.Surv_PB(:,2),'ok');
    hold on;
    plot(handles.Hist.Surv_PB(:,1),handles.Hist.Surv_PB(:,3),'g', 'LineWidth', 1);
    plot(handles.Hist.Surv_PB(:,1),handles.Hist.Surv_PB(:,4),'r', 'LineWidth', 1);
    plot(handles.Hist.Surv_PB(:,1),handles.Hist.Surv_PB(:,5),'b', 'LineWidth', 1);
    plot(handles.Hist.Surv_PB(:,1),handles.Hist.Surv_PB(:,6),'b--', 'LineWidth', 1);
    hold off;
    xlabel('Time [s]');
    ylabel('Counts');
    title({['Fit of the survival probability',...
        ' - Measured [partial] C_{eq} = ' num2str(handles.FitPar.Surv_PB(1,7),3), ' \pm ',...
        num2str(handles.FitPar.Surv_PB(2,7),2)],...
        ['One Component fit k_{off} = ',num2str(handles.FitPar.Surv_PB(1,1),3),...
        ' \pm ', num2str(handles.FitPar.Surv_PB(2,7),2), ' s^{-1}; C_{eq} = ', num2str(handles.FitPar.Surv_PB(1,2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(handles.FitPar.Surv_PB(1,3),3),...
        ' \pm ', num2str(handles.FitPar.Surv_PB(2,3),2), ' s^{-1}'], ...
        ['k_{off, 2} = ', num2str(handles.FitPar.Surv_PB(1,4),3),...
        ' \pm ', num2str(handles.FitPar.Surv_PB(2,4),2), ' s^{-1}; f_1 = ', num2str(handles.FitPar.Surv_PB(1,5),3),...
        '\pm', num2str(handles.FitPar.Surv_PB(2,5),2), 'C_{eq} = ',  num2str(handles.FitPar.Surv_PB(1,6),2)]});
    legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit');
    box on;
    set(handles.PieChart,'Enable','on');
    
    %Display Parameters
    set(handles.frame2frameMaxText,'String',['Max. Frame to Frame Jump: ', num2str(handles.Parameters.ThreshL), ' um']);
    set(handles.End2EndMaxText,'String',['Max. End to End Jump: ', num2str(handles.Parameters.ThreshH), ' um']);
    set(handles.NminText,'String',['Min. Bound Frames (Nmin): ', num2str(handles.Parameters.minBoundFrames)]);
    guidata(hObject,handles);
end


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function PieChart_Callback(hObject, eventdata, handles)
% hObject    handle to PieChart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fractions,resTimes] = PieChartsFromResTimeAnalysis(handles);

frac = fractions';
res = resTimes';

data = [frac(:,1:2),res(:,1)];
if size(frac,2) == 3
    data = [data,frac(:,3),res(:,2)];
end

num2clip(data);
msgbox('Fractions and Residence Times have been copied to the clipboard');


% --------------------------------------------------------------------
function MergeAnalyze_Callback(hObject, eventdata, handles)
% hObject    handle to MergeAnalyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileNames,PathName,FilterIndex] = uigetfile('.mat',...
    'Select Mat files for the analysis of the residence time',...
    'MultiSelect','on');
if ~iscell(FileNames)
    if FileNames == 0
        FileNames = [];
    else
         FileNames = {FileNames};
    end
   
end

if ~isempty(FileNames) 
    if isfield(handles,'Parameters')
        if isfield(handles.Parameters, 'Used')
            if isfield(handles.Parameters.Used, 'Bound');
                
                defPar(1) = handles.Parameters.Used.Bound.ThreshL;
                defPar(2) = handles.Parameters.Used.Bound.ThreshH;
                defPar(3) = handles.Parameters.Used.Bound.minBoundFrames;
                defPar(4) = handles.Parameters.Used.Acquisition.frameTime;
                
            else
                defPar(1) = handles.Parameters.Bound.ThreshL;
                defPar(2) = handles.Parameters.Bound.ThreshH;
                defPar(3) = handles.Parameters.Bound.minBoundFrames;
                defPar(4) = handles.Parameters.Acquisition.frameTime;
            end
        else
            defPar(1) = 0.22;
            defPar(2) = 0.315;
            defPar(3) = 2;
            defPar(4) = 0.2;
            
        end
    else
        defPar(1) = 0.22;
        defPar(2) = 0.315;
        defPar(3) = 2;
        defPar(4) = 0.2;
    end
        
    
    % Ask if survival distr or residence time hist needs to be analyzed;
    answer = questdlg({'Do you want to analyze';...
        'the residence time histogram'; ...
        'or the survival time distribution?'},'Analyze?', ...
        'Residence Time Histogram', 'Survival Time Distribution',...
        'Residence Time Histogram');
    
    switch answer;
        case 'Residence Time Histogram'
            AnalysisFlag = 1;
            
            
            
        case 'Survival Time Distribution'
            AnalysisFlag = 2;
            
    end
    
    % Get a dialog input to set the parameters;
    %     if AnalysisFlag == 1
    %         defaults = {num2str(defPar(1)), num2str(defPar(2)), ...
    %             num2str(defPar(3)), '1', num2str(defPar(4))};
    %         prompt = {'Maximum Jump between consecutive frames',...
    %             'Maximum end to end distance',...
    %             'Minimum length of bound tracks',...
    %             'Histogram binning(1 = no binning)',...
    %             'Frame Time'};
    %         dlgtitle = 'Parameters for the residence time histogram';
    %         AnalysisParam = inputdlg(prompt,dlgtitle, 1, defaults);
    %
    %     elseif AnalysisFlag == 2
    defaults = {num2str(defPar(1)), num2str(defPar(2)), ...
        num2str(defPar(3)), '0', num2str(defPar(4))};
    prompt = {'Maximum Jump between consecutive frames',...
        'Maximum end to end distance',...
        'Minimum length of bound tracks',...
        'Points in the histogram after log sampling (0 = No Log Sampling)',...
        'Frame Time'};
    dlgtitle = 'Parameters for the Survival Time distribution';
    AnalysisParam = inputdlg(prompt,dlgtitle, 1, defaults);
    
    %     end
    
    % Send parameters and file names to the residence time GUI
    handles.Flag = AnalysisFlag;
    handles.PathName = PathName;
    handles.FileNames = FileNames;
    AnalysisS = AnalysisParam;

    handles.Parameters = [];
    handles.Parameters.ThreshL = str2num(AnalysisS{1});
    handles.Parameters.ThreshH = str2num(AnalysisS{2});
    handles.Parameters.minBoundFrames = str2num(AnalysisS{3});
    handles.Parameters.bin = str2num(AnalysisS{4});
    handles.FrameTime = str2num(AnalysisS{5});




    % If only one file has been selected change the fileNames variable
    % to a cell array

    if ~iscell(handles.FileNames)
        handles.FileNames = {handles.FileNames};
    end;
    nFiles = length(handles.FileNames);

    % get the tracks of the selected files
    % open each mat-files and retrieve useful information
    ROIChoice = {};
    ROIorClass = -1; % -1: not chosen, 0: ROI name, 1: Class Name
    bgVec = [];
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

        if isfield(Temp.Results, 'PreAnalysis') && ...
                isfield(Temp.Results.PreAnalysis, 'Tracks_um');
            Tracks{i} = Temp.Results.PreAnalysis.Tracks_um;
            NParticles{i} = Temp.Results.PreAnalysis.NParticles;
            if isfield(Temp.Results,'Process') %check that the strcutures exist
                if isfield(Temp.Results.Process,'ROIlabel')
                    if isfield(Temp.Results.Process,'AllROIClasses') && ROIorClass == -1
                        ROIClass =  questdlg('Do you want to separate data based on ROI names or Class Names?','ROI or Class','ROI','Class','Class');
                        if strcmp(ROIClass,'ROI');
                            ROIorClass = 0;
                        else
                            ROIorClass = 1;
                        end
                    elseif ~isfield(Temp.Results.Process,'AllROIClasses')
                        ROIorClass = 0;
                    end
                    if ROIorClass == 0
                        roiLabels = Temp.Results.Process.ROIlabel;
                    elseif ROIorClass == 1
                        if ~isfield(Temp.Results.Process,'AllROIClasses')
                            errordlg(['File does not contain Class data: ,' handles.FileNames{i}]);
                        else
                            roiLabels = Temp.Results.Process.AllROIClasses;
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
                                if size(Tracks{i},2) < 8
                                    Tracks{i} = Tracks{i}(Tracks{i}(:,5) == ROIidx,:);
                                else
                                    Tracks{i} = Tracks{i}(Tracks{i}(:,9) == ROIidx,:);
                                end
                                NParticles{i} = [NParticles{i}(:,1) NParticles{i}(:,ROIidx+1)];
                                ROIChoice{end+1,1} = Temp.Results.Process.ROIlabel{ROIidx};
                            else
                                roiAct = [];
                                tracks_tmp = Tracks{i};
                                Tracks{i} = [];
                                for m = 1:size(Temp.Results.Process.ROIClass);
                                    if strcmpi(roiLabels{ROIidx,:},Temp.Results.Process.ROIClass{m,:})
                                        if size(Tracks{i},2) < 8
                                            Tracks{i} = [Tracks{i}; tracks_tmp(tracks_tmp(:,5) == m,:)];
                                        else
                                            Tracks{i} = [Tracks{i}; tracks_tmp(tracks_tmp(:,9) == m,:)];
                                        end
                                        roiAct = [roiAct;m];
                                    end

                                end
                                tmp = NParticles{i}(:,roiAct+1);

                                NParticles{i} = [NParticles{i}(:,1) sum(tmp,2)];
                                ROIChoice{end+1,1} = Temp.Results.Process.AllROIClasses{ROIidx};
                            end

                        else
                            tmp = NParticles{i}(:,2:size(Temp.Results.Process.ROIlabel,1)+1);
                            NParticles{i} = [NParticles{i}(:,1) sum(tmp,2)];
                            ROIChoice{end+1,1} = 'All';
                        end
                    else
                        ROIChoice{1} = roiLabels{1};
                    end
                else
                    ROIChoice{1} = '';
                end
            else
                ROIChoice{1} = '';
            end


            if ~isempty(Tracks{i})
                if size(Tracks{i},2) < 8
                    [ImmTracks, Dummy] = calculateImmobileTracks...
                        (Tracks{i}, handles.Parameters.ThreshL,...
                        handles.Parameters.minBoundFrames, handles.Parameters.ThreshH, 0);


                    if isempty(ImmTracks)
                        TrackLengthHist{i} = [0 0];
                    else

                        % Identify bound molecules for each of the movies;

                        % Calculate survival histogram;
                        TrackLengthHist{i} = ...
                            calculateTrackLength(ImmTracks, ...
                            handles.FrameTime,handles.Parameters.minBoundFrames);



                    end
                else
                    ImmTracks = Tracks{i}(Tracks{i}(:,2) >= handles.Parameters.minBoundFrames,:);
                    ImmTracks = ImmTracks(ImmTracks(:,3) <= handles.Parameters.ThreshL,:);
                    if ~isempty(ImmTracks)
                        TrackLength = ImmTracks(:,2);
                        LongestTrack = max(TrackLength);
                        ShortestTrack = handles.Parameters.minBoundFrames;

                        TrackLengthHist{i} = zeros(LongestTrack-ShortestTrack + 1,2);
                        TrackLengthHist{i}(:,1) = (ShortestTrack-1: LongestTrack)* handles.FrameTime;
                        %         TrackLengthHist(:,2) = hist(TrackLength,(ShortestTrack:LongestTrack));
                        for  j = 1:length(TrackLength);
                            TrackLengthHist{i}(1:TrackLength(j)- ShortestTrack + 1,2) = ...
                                TrackLengthHist{i}(1:TrackLength(j) - ShortestTrack + 1,2)+1;
                        end;
                    end
                end
            end


        else
            errordlg(['Some of the files have not been preprocessed: ',handles.PathName,handles.FileNames{i}]);
            return
        end
    end
    %Calculate the avg image intensity over time
    % avgInt = sum(pxIntesity)./sum(npx);

    % Calculate the total number of particles;

    % find the data containing the lowest number of frames;
    TimePoints = [];


    nMax = length(NParticles{1}(:,1));
    TimePoints = NParticles{1}(:,1);
    for i = 1:nFiles

        if length(NParticles{i}(:,1)) < nMax;

            TimePoints = NParticles{i}(:,1);
            nMax = length(TimePoints);
        end
    end

    if isempty(TimePoints)
        errordlg('Too few bound particles to produce an histogram');
        return
    end

    % Accumulate the histogram
    Hist_Matrix = zeros(nMax,nFiles);

    for i = 1:nFiles
        Hist_Matrix(1:nMax,i) = NParticles{i}(1:nMax,2);
    end

    CumNParticles(:,1) = (TimePoints-1) * handles.FrameTime;
    CumNParticles(:,2) = sum(Hist_Matrix, 2);

    %Calculate the avg intensity


    % Accumulate the different survival histograms;

    % find the histogram containing the longest track
    nMax = 1;
    TimePoints = [];
    totalTracks = 0;
    for i = 1:nFiles
        trInd = unique(Tracks{i}(:,4));
        totalTracks = totalTracks + length(trInd);

        if length(TrackLengthHist{i}(:,1)) > nMax;

            TimePoints = TrackLengthHist{i}(:,1);
            nMax = length(TimePoints);
        end
    end
    handles.totalTracks = totalTracks;
    if isempty(TimePoints)
        errordlg('Too few bound particles to produce an histogram');
        return
    end


    % Accumulate the histogram
    Hist_Matrix = zeros(nMax,nFiles);

    for i = 1:nFiles
        lengthHist = length(TrackLengthHist{i}(:,1));
        Hist_Matrix(1:lengthHist,i) = TrackLengthHist{i}(:,2);
    end

    CumTrackLengthHist(:,1) = TimePoints;
    CumTrackLengthHist(:,2) = sum(Hist_Matrix, 2);

    % Copy Histogram to handles
    handles.CumHist = CumTrackLengthHist;


    % Calculate the fraction of bound molecules.

    % First calculate the residence time histogram
    ResTimeHist(:,1)= CumTrackLengthHist(:,1);
    ResTimeHist(:,2) = ...
        [(CumTrackLengthHist(1:end-1,2) - CumTrackLengthHist(2:end,2));CumTrackLengthHist(end,2)];

    % Then calculate the total number of bound spots
    TotalBoundMolecules = sum(ResTimeHist(:,2).*ResTimeHist(:,1)/handles.FrameTime);
    handles.BoundMolecules = TotalBoundMolecules;
    % and the total number of spots
    % TotalMolecules = sum(CumNParticles(:,2));
    TotalMolecules = 0;
    for i = 1:nFiles
        TotalMolecules = TotalMolecules + sum(NParticles{i}(:,2));
    end
    handles.TotalMolecules = TotalMolecules;
    % Finally divide the two.
    PartialBoundFraction = TotalBoundMolecules/TotalMolecules;

    BFerror = (sqrt(TotalBoundMolecules)/TotalBoundMolecules + ...
        sqrt(TotalMolecules)/TotalMolecules)*PartialBoundFraction;


    % NOTE: THIS IS ONLY THE PARTIAL BOUND FRACTION
    % BECAUSE WE ARE LOOKING ONLY AT PARTICLES BOUND
    % FOR MORE THAN Nmin. We HAVE TO FIT THE HISTOGRAM OF RESIDENCE
    % TIMES AND EXTRAPOLATE THE TRUE VALUE OF THE BOUND FRACTION

    % Photobleaching correction




    % Fit a double exponential to the photobleaching curve
    disp(' ')
    disp('____________________________________')
    disp('Estimating bleaching characteristics')
    disp('____________________________________')

    [BleachRates, Dummy, CumNParticles(:,3:6)] =...
        ExpDecay_2Cmp_fit(CumNParticles, [1 0.1]);
    disp(['Bleach Rate 1: ', num2str(BleachRates(1), 3), ' s^-1'])
    disp(['Bleach Rate 2: ', num2str(BleachRates(2), 3), ' s^-1'])
    disp(['Fraction 1: ', num2str(BleachRates(3), 3)])
    disp('____________________________________')
    disp(' ')

    CumNParticles = CumNParticles(:,1:4);
    % Plot the photobleaching curve
%     PhotobAxes = findobj ('Tag', 'photob_axes');
    axes(handles.photob_axes);






    plot(CumNParticles(:,1), CumNParticles(:,2)/max(CumNParticles(:,4)), '+', 'MarkerEdgeColor', [0.5 0.5 0.5]');
    hold on
    plot(CumTrackLengthHist(:,1), CumTrackLengthHist(:,2)/max(CumTrackLengthHist(:,2)),'ok');
    plot(CumNParticles(:,3), CumNParticles(:,4)/max(CumNParticles(:,4)),'r');
    hold off;

    box on;

    title({'Photobleaching kinetics:', ['k_{b1} = ', num2str(BleachRates(1), 3),...
        ' s^{-1}, k_{b2} = ', num2str(BleachRates(2), 3), ' s^{-1}, f_1 = ',  ...
        num2str(BleachRates(3), 3)]})

    hold off;
    legend('Photobleaching Decay', 'Bound molecules Decay');

    xlabel('Time [s]');
    ylabel('Normalized Counts');
    % Copy CumNparticles to handles
    handles.NParticles = CumNParticles;
    handles.BleachRates = BleachRates;

    CumTrackLengthHist_PB = CumTrackLengthHist;
    CumTrackLengthHist_PB(:,2) = CumTrackLengthHist_PB(:,2)./ ...
        (BleachRates(3)*exp(-BleachRates(1).* CumTrackLengthHist_PB(:,1)) + ...
        (1-BleachRates(3))*exp(-BleachRates(2).* CumTrackLengthHist_PB(:,1)));
    ResTimeHist_PB(:,1)= CumTrackLengthHist_PB(:,1);
    ResTimeHist_PB(:,2) = ...
        [(CumTrackLengthHist_PB(1:end-1,2) - CumTrackLengthHist_PB(2:end,2));CumTrackLengthHist_PB(end,2)];

    %only use values where there are particles with that residence time for curve-fitting
    ind2fit = find(ResTimeHist(:,2) > 0);


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
    deltaT = ResTimeHist(2,1) - ResTimeHist(1,1);
    % Normalize the Survival probability for the partial bound fraction
    % if size(Tracks{1},2) < 8
    CumTrackLengthHist(:,2) = CumTrackLengthHist(:,2)/CumTrackLengthHist(1,2)...
        .*PartialBoundFraction;

    CumTrackLengthHist_PB(:,2) = CumTrackLengthHist_PB(:,2)/CumTrackLengthHist_PB(1,2)...
        .*PartialBoundFraction;
    % Normalize the residence time histogram for the partial bound fraction
    ResTimeHist(:,2) = ResTimeHist(:,2)/(deltaT*sum(ResTimeHist(:,2)))*...
        PartialBoundFraction;

    ResTimeHist_PB(:,2) = ResTimeHist_PB(:,2)/(deltaT*sum(ResTimeHist_PB(:,2)))*...
        PartialBoundFraction;
    % end


    % handles.Hist.Surv(:,1:2) = CumTrackLengthHist;
    % 
    % handles.Hist.Surv_PB(:,1:2) = CumTrackLengthHist_PB;
    % 
    % handles.Hist.Res(:,1:2) = ResTimeHist;
    % 
    % handles.Hist.Res_PB(:,1:2) = ResTimeHist_PB;
    handles.PartBoundFrac = PartialBoundFraction;
    handles.BFerror = BFerror;

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
            binFactor = floor(length(ResTimeHist(:,1))/nBins);
    %     binFactor = handles.Parameters.bin + 1;
    %     nBins = floor(length(ResTimeHist(:,1))/binFactor);
    %     if strcmp(answer,'Yes')
            TimeBins = reshape(ResTimeHist(1:nBins*binFactor,1), binFactor, nBins);
            CountsBins = reshape(ResTimeHist(1:nBins*binFactor,2), binFactor, nBins);
    %     else
            TimeBins_PB = reshape(ResTimeHist_PB(1:nBins*binFactor,1), binFactor, nBins);
            CountsBins_PB = reshape(ResTimeHist_PB(1:nBins*binFactor,2), binFactor, nBins);
    %     end


            resTimeHist_Binned = [];
            resTimeHist_Binned(:,1) = mean(TimeBins);
            resTimeHist_Binned(:,2) = sum(CountsBins);
            resTimeHist_Binned_PB = [];
            resTimeHist_Binned_PB(:,1) = mean(TimeBins_PB);
            resTimeHist_Binned_PB(:,2) = sum(CountsBins_PB);
        else
            resTimeHist_Binned = ResTimeHist;
            resTimeHist_Binned_PB = ResTimeHist_PB;
        end

        % Remove values < 0 from the histogram for fitting only
        idx = find(resTimeHist_Binned(:,2) < 0);
        resTimeHist_Binned(idx,2) = 0;
        idx = find(resTimeHist_Binned(:,2) <= 0);
        resTimeHist_Binned2 = resTimeHist_Binned;
        resTimeHist_Binned2(idx,:) = [];
        % Fit the binned residence time histogram with exponentials;
        [fitpar,espSigma, ~]= ExpDecay_fit_resTime(resTimeHist_Binned2, 1 );
        fit = [resTimeHist_Binned(:,1),ExpDecay_fun_resTime(fitpar,resTimeHist_Binned(:,1))];
        [fitpar2,espSigma2, ~]= ExpDecay_2Cmp_fit_resTime(resTimeHist_Binned2, [0.1 0.01]);

        Esp_Coef_1 = [fitpar2(1) fitpar2(4)*fitpar2(3)];
        Esp_Fit(:,1) = ExpDecay_fun_resTime(Esp_Coef_1,resTimeHist_Binned(:,1));

        Esp_Coef_2 = [fitpar2(2) fitpar2(4)*(1 - fitpar2(3))];
        Esp_Fit(:,2) = ExpDecay_fun_resTime(Esp_Coef_2,resTimeHist_Binned(:,1));

        fit2 = [resTimeHist_Binned(:,1),ExpDecay_2Cmp_fun_resTime(fitpar2,resTimeHist_Binned(:,1)),Esp_Fit];

        handles.Hist.Res(:, 1:2) = resTimeHist_Binned;
        handles.Hist.Res(:,3) = fit(:,2);
        handles.Hist.Res(:,4) = fit2(:,2);
        handles.Hist.Res(:,5) = fit2(:,3);
        handles.Hist.Res(:,6) = fit2(:,4);

        [~, pvalRes,FstatRes] = FtestModelCompare(resTimeHist_Binned(:,2),fit(:,2),fit2(:,2),2,4);

        handles.FitPar.Res = [fitpar, fitpar2, PartialBoundFraction, pvalRes; espSigma, espSigma2, BFerror,FstatRes];
        idx = find(resTimeHist_Binned_PB(:,2) < 0);
        resTimeHist_Binned_PB(idx,2) = 0;
        idx = find(resTimeHist_Binned_PB(:,2) <= 0);
        resTimeHist_Binned_PB2 = resTimeHist_Binned_PB;
        resTimeHist_Binned_PB2(idx,:) = [];
        % Fit the binned residence time histogram with exponentials;
        [fitpar_PB,espSigma_PB, fit]= ExpDecay_fit_resTime(resTimeHist_Binned_PB2, 1 );
        fit_PB = [resTimeHist_Binned(:,1),ExpDecay_fun_resTime(fitpar_PB,resTimeHist_Binned_PB(:,1))];

        [fitpar2_PB,espSigma2_PB, fit2]= ExpDecay_2Cmp_fit_resTime(resTimeHist_Binned_PB2, [0.1 0.01]);

        Esp_Coef_1 = [fitpar2_PB(1) fitpar2_PB(4)*fitpar2_PB(3)];
        Esp_Fit(:,1) = ExpDecay_fun_resTime(Esp_Coef_1,resTimeHist_Binned_PB(:,1));

        Esp_Coef_2 = [fitpar2_PB(2) fitpar2_PB(4)*(1 - fitpar2_PB(3))];
        Esp_Fit(:,2) = ExpDecay_fun_resTime(Esp_Coef_2,resTimeHist_Binned_PB(:,1));

        fit2_PB = [resTimeHist_Binned_PB(:,1),ExpDecay_2Cmp_fun_resTime(fitpar2_PB,resTimeHist_Binned_PB(:,1)),Esp_Fit];
        handles.Hist.Res_PB(:, 1:2) = resTimeHist_Binned_PB;
        handles.Hist.Res_PB(:,3) = fit_PB(:,2);
        handles.Hist.Res_PB(:,4) = fit2_PB(:,2);
        handles.Hist.Res_PB(:,5) = fit2_PB(:,3);
        handles.Hist.Res_PB(:,6) = fit2_PB(:,4);
        [~, pvalRes_PB,FstatRes_PB] = FtestModelCompare(resTimeHist_Binned_PB(:,2),fit_PB(:,2),fit2_PB(:,2),2,4);
        handles.FitPar.Res_PB = [fitpar_PB, fitpar2_PB, PartialBoundFraction, pvalRes_PB; espSigma_PB, espSigma2_PB, BFerror, FstatRes_PB];


        HistAxes = findobj ('Tag', 'hist_axes');
        axes(HistAxes);
        if strcmp(answer,'Yes');
            resTimeHist_Binned = resTimeHist_Binned_PB;
            fit = fit_PB;
            fit2 = fit2_PB;
            fitpar = fitpar_PB;
            fitpar2 = fitpar2_PB;
            espSigma = espSigma_PB;
            espSigma2 = espSigma2_PB;
            pvalRes = pvalRes_PB;
            FstatRes = FstatRes_PB;
        end
        if handles.Flag == 1
            bar(resTimeHist_Binned(:,1), resTimeHist_Binned(:,2),'EdgeColor', [0 0 0], ...
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
                ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction,3), ' \pm ',...
                num2str(BFerror,2)],...
                ['One Component fit k_{off} = ',num2str(fitpar(1),3),...
                ' \pm ', num2str(espSigma(1),2), ' s^{-1}; C_{eq} = ', num2str(fitpar(2),2)],...
                ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
                ' \pm ', num2str(espSigma2(1),2), ' s^{-1}'], ...
                ['k_{off, 2} = ', num2str(fitpar2(2),3),...
                ' \pm ', num2str(espSigma2(2),2), ' s^{-1}; f_1 = ', num2str(fitpar2(3),3),...
                '\pm', num2str(espSigma2(3),2), ' C_{eq} = ',  num2str(fitpar2(4),2)]});
            legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit','Component 1', 'Component 2');
            box on;
            set(handles.FstatText,'String',['F stat: ', num2str(FstatRes,4)]);
            set(handles.pValText,'String',['p-value: ', num2str(pvalRes,4)]);

        end




    % end

    % if handles.Flag == 2;
    %     set(handles.What2Plot,'Value',2);
        % Logarithmic sampling of the survival prob
        if handles.Parameters.bin ~= 0

    %         if strcmp(answer,'Yes')
                TimePoints_PB = logspace(log10(min(CumTrackLengthHist_PB(:,1))),...
                    log10(max(CumTrackLengthHist_PB(:,1))),handles.Parameters.bin);
    %         else
                TimePoints = logspace(log10(min(CumTrackLengthHist(:,1))),...
                    log10(max(CumTrackLengthHist(:,1))),handles.Parameters.bin);
    %         end
            TrLog =[];

            TimePoints = round(TimePoints*1000)/1000;
            TimePoints_PB = round(TimePoints_PB*1000)/1000;
            for i = 1:length(TimePoints)

                idx = find(CumTrackLengthHist(:,1) <= TimePoints(i)+0.0001,1,'last');
                TrLog(i,:) = CumTrackLengthHist(idx,:);
            end
            for i = 1:length(TimePoints_PB)

                idx = find(CumTrackLengthHist_PB(:,1) <= TimePoints_PB(i)+0.0001,1,'last');
                TrLog_PB(i,:) = CumTrackLengthHist_PB(idx,:);
            end

        else
            TrLog = CumTrackLengthHist;
            TrLog_PB = CumTrackLengthHist_PB;
        end

        tpoint1 = TrLog(1,1);
        TrLog(:,1) = TrLog(:,1) - tpoint1;


    %     ind2fit = 1:length(TrLog);
        % Fit exponentials to the Survival probability
        ind2fit = 1:size(TrLog,1);
        if handles.Parameters.bin == 0

            [fitpar,espSigma, ~]= ExpDecay_fit(TrLog(ind2fit,:), 1/(tpoint1*10));
            fit = [TrLog(:,1),ExpDecay_fun(fitpar,TrLog(:,1))];
            [fitpar2,espSigma2, ~]= ExpDecay_2Cmp_fit(TrLog(ind2fit,:), [1/(tpoint1*10) 1/(tpoint1*100)]);
            Esp_Coef_1 = [fitpar2(1) fitpar2(4)*fitpar2(3)];
            Esp_Fit(:,1) = ExpDecay_fun(Esp_Coef_1,TrLog(:,1));

            Esp_Coef_2 = [fitpar2(2) fitpar2(4)*(1 - fitpar2(3))];
            Esp_Fit(:,2) = ExpDecay_fun(Esp_Coef_2,TrLog(:,1));
            fit2 = [TrLog(:,1),ExpDecay_2Cmp_fun(fitpar2,TrLog(:,1)),Esp_Fit];
        else
            [fitpar,espSigma, fit]= ExpDecay_fit(TrLog(ind2fit,:), 1/(tpoint1*10));
            [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit(TrLog(ind2fit,:), [1/(tpoint1*10) 1/(tpoint1*100)]);
        end

        fit(:,1) = fit(:,1) + tpoint1;
        fit2(:,1) = fit2(:,1) + tpoint1;
        TrLog(:,1) = TrLog(:,1) + tpoint1;

        tpoint1 = TrLog_PB(1,1);
        TrLog_PB(:,1) = TrLog_PB(:,1) - tpoint1;



        % Fit exponentials to the Survival probability
        if handles.Parameters.bin == 0
            [fitpar_PB,espSigma_PB, ~]= ExpDecay_fit(TrLog_PB(ind2fit,:), 1/(tpoint1*10));
            fit_PB = [TrLog_PB(:,1),ExpDecay_fun(fitpar_PB,TrLog_PB(:,1))];
            [fitpar2_PB,espSigma2_PB, ~]= ExpDecay_2Cmp_fit(TrLog_PB(ind2fit,:), [1/(tpoint1*10) 1/(tpoint1*100)]);
            Esp_Coef_1 = [fitpar2_PB(1) fitpar2_PB(4)*fitpar2_PB(3)];
            Esp_Fit(:,1) = ExpDecay_fun(Esp_Coef_1,TrLog_PB(:,1));

            Esp_Coef_2 = [fitpar2_PB(2) fitpar2_PB(4)*(1 - fitpar2_PB(3))];
            Esp_Fit(:,2) = ExpDecay_fun(Esp_Coef_2,TrLog_PB(:,1));
            fit2_PB = [TrLog_PB(:,1),ExpDecay_2Cmp_fun(fitpar2_PB,TrLog_PB(:,1)),Esp_Fit];
        else
            [fitpar_PB,espSigma_PB, fit_PB]= ExpDecay_fit(TrLog_PB(ind2fit,:), 1/(tpoint1*10));
            [fitpar2_PB,espSigma2_PB, fit2_PB]= ExpDecay_2Cmp_fit(TrLog_PB(ind2fit,:), [1/(tpoint1*10) 1/(tpoint1*100)]);
        end
        fit_PB(:,1) = fit_PB(:,1) + tpoint1;
        fit2_PB(:,1) = fit2_PB(:,1) + tpoint1;
        TrLog_PB(:,1) = TrLog_PB(:,1) + tpoint1;

        handles.Hist.Surv(:, 1:2) = TrLog;
        handles.Hist.Surv(:,3) = fit(:,2);
        handles.Hist.Surv(:,4) = fit2(:,2);
        handles.Hist.Surv(:,5) = fit2(:,3);
        handles.Hist.Surv(:,6) = fit2(:,4);
        [~, pvalSurv, FstatSurv] = FtestModelCompare(TrLog(:,2),fit(:,2),fit2(:,2),2,4);
        handles.FitPar.Surv = [fitpar, fitpar2, PartialBoundFraction, pvalSurv; espSigma, espSigma2, BFerror, FstatSurv];

        handles.Hist.Surv_PB(:, 1:2) = TrLog_PB;
        handles.Hist.Surv_PB(:,3) = fit_PB(:,2);
        handles.Hist.Surv_PB(:,4) = fit2_PB(:,2);
        handles.Hist.Surv_PB(:,5) = fit2_PB(:,3);
        handles.Hist.Surv_PB(:,6) = fit2_PB(:,4);
        [~, pvalSurv_PB, FstatSurv_PB] = FtestModelCompare(TrLog_PB(:,2),fit_PB(:,2),fit2_PB(:,2),2,4);
        handles.FitPar.Surv_PB = [fitpar_PB, fitpar2_PB, PartialBoundFraction, pvalSurv_PB; espSigma_PB, espSigma2_PB, BFerror, FstatSurv_PB];
        if strcmp(answer,'Yes');
            TrLog = TrLog_PB;
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
            HistAxes = findobj ('Tag', 'hist_axes');
            axes(HistAxes);
            plot(TrLog(:,1), TrLog(:,2),'ok');
            hold on;
            plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
            plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
            plot(fit2(:,1),fit2(:,3),'b', 'LineWidth', 1);
            plot(fit2(:,1),fit2(:,4),'b--', 'LineWidth', 1);
            hold off;
            xlabel('Time [s]');
            ylabel('Counts');
            title({['Fit of the survival probability',...
                ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction,3), ' \pm ',...
                num2str(BFerror,2)],...
                ['One Component fit k_{off} = ',num2str(fitpar(1),3),...
                ' \pm ', num2str(espSigma(1),2), ' s^{-1}; C_{eq} = ', num2str(fitpar(2),2)],...
                ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
                ' \pm ', num2str(espSigma2(1),2), ' s^{-1}'], ...
                ['k_{off, 2} = ', num2str(fitpar2(2),3),...
                ' \pm ', num2str(espSigma2(2),2), ' s^{-1}; f_1 = ', num2str(fitpar2(3),3),...
                '\pm', num2str(espSigma2(3),2), 'C_{eq} = ',  num2str(fitpar2(4),2)]});
            legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit');
            box on;
            set(handles.FstatText,'String',['F stat: ', num2str(FstatSurv,4)]);
            set(handles.pValText,'String',['p-value: ', num2str(pvalSurv,4)]);
        end





    set(handles.totalTracksText,'String',['Total Tracks: ',num2str(handles.totalTracks)]);
    set(handles.boundParticlesText,'String',['Bound Particles: ',num2str(handles.BoundMolecules)]);
    set(handles.totalParticlesText,'String',['Total Particles: ',num2str(handles.TotalMolecules)]);

    figname = get(gcf,'Name');
    figname = [figname ' ' ROIChoice{1}];
    set(gcf,'Name',figname);
    set(handles.PieChart,'Enable','on');

end

    % Update handles structure
    guidata(hObject, handles);


% --------------------------------------------------------------------
function comparePar_Callback(hObject, eventdata, handles)
% hObject    handle to comparePar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CompareParCI_GUI();
