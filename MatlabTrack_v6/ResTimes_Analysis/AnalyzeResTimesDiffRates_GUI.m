function varargout = AnalyzeResTimesDiffRates_GUI(varargin)
% ANALYZERESTIMESDIFFRATES_GUI MATLAB code for AnalyzeResTimesDiffRates_GUI.fig
%      ANALYZERESTIMESDIFFRATES_GUI, by itself, creates a new ANALYZERESTIMESDIFFRATES_GUI or raises the existing
%      singleton*.
%
%      H = ANALYZERESTIMESDIFFRATES_GUI returns the handle to a new ANALYZERESTIMESDIFFRATES_GUI or the handle to
%      the existing singleton*.
%
%      ANALYZERESTIMESDIFFRATES_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYZERESTIMESDIFFRATES_GUI.M with the given input arguments.
%
%      ANALYZERESTIMESDIFFRATES_GUI('Property','Value',...) creates a new ANALYZERESTIMESDIFFRATES_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AnalyzeResTimesDiffRates_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AnalyzeResTimesDiffRates_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AnalyzeResTimesDiffRates_GUI

% Last Modified by GUIDE v2.5 16-May-2017 10:25:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @AnalyzeResTimesDiffRates_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @AnalyzeResTimesDiffRates_GUI_OutputFcn, ...
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


% --- Executes just before AnalyzeResTimesDiffRates_GUI is made visible.
function AnalyzeResTimesDiffRates_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AnalyzeResTimesDiffRates_GUI (see VARARGIN)

% Choose default command line output for AnalyzeResTimesDiffRates_GUI
handles.output = hObject;

% Import Parameters.
handles.Flag = varargin{1};
handles.PathName1 = varargin{2};
handles.FileNames1 = varargin{3};
handles.PathName2 = varargin{4};
handles.FileNames2 = varargin{5};
AnalysisS = varargin{6};

handles.Parameters = [];
handles.Parameters.ThreshL = str2double(AnalysisS{1});
handles.Parameters.ThreshH = str2double(AnalysisS{2});
handles.Parameters.minBoundFrames = str2double(AnalysisS{3});
handles.Parameters.bin = str2double(AnalysisS{4});
handles.FrameTime1 = str2double(AnalysisS{5});
handles.FrameTime2 = str2double(AnalysisS{6});




% If only one file has been selected change the fileNames variable
% to a cell array

if ~iscell(handles.FileNames1)
    handles.FileNames1 = {handles.FileNames1};
end;
nFiles1 = length(handles.FileNames1);
if ~iscell(handles.FileNames2)
    handles.FileNames2 = {handles.FileNames2};
end;
nFiles2 = length(handles.FileNames2);

% get the tracks of the selected files
% open each mat-files and retrieve useful information
ROIChoice = {};
bgVec = [];
for i = 1:nFiles1
    
    ROIidx = 0;
    Temp = load([handles.PathName1,handles.FileNames1{i}],'Results');
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
                if size(Temp.Results.Process.ROIlabel,1) > 1
                    if ~isempty(ROIChoice) %if a named ROI was selected before check to see if an ROI with that name exists in this dataset
                        if ~strcmp(ROIChoice{1},'All')
                            for j = 1:length(ROIChoice)
                                for k = 1:length(Temp.Results.Process.ROIlabel)
                                    if strcmpi(ROIChoice{j},Temp.Results.Process.ROIlabel{k})
                                        ROIidx = k;
                                        break
                                    end
                                end
                            end
                        else
                            ROIidx = size(Temp.Results.Process.ROIlabel,1) + 1; % If all was selected before, just set the index to something greater than the number of ROIs
                        end
                    end
                    
                    if ROIidx == 0 %if nothing has been found automatically, provide the user with a dialog box to select the ROI
                        
                        ROIstring = Temp.Results.Process.ROIlabel;
                        ROIstring{end+1,1} = 'All';
                        ROIidx = ROIchooseDlg(ROIstring);
                        
                    end
                    %Update the Tracks & NParticles data
                    if ROIidx <= size(Temp.Results.Process.ROIlabel,1)
                        if size(Tracks{i},2) < 8
                            Tracks{i} = Tracks{i}(Tracks{i}(:,5) == ROIidx,:);
                        else
                            Tracks{i} = Tracks{i}(Tracks{i}(:,9) == ROIidx,:);
                        end
                        NParticles{i} = [NParticles{i}(:,1) NParticles{i}(:,ROIidx+1)];
                        ROIChoice{end+1,1} = Temp.Results.Process.ROIlabel{ROIidx};
                    else
                        tmp = NParticles{i}(:,2:size(Temp.Results.Process.ROIlabel,1)+1);
                        NParticles{i} = [NParticles{i}(:,1) sum(tmp,2)];
                        ROIChoice{end+1,1} = 'All';
                    end
                else
                    ROIChoice{1} = Temp.Results.Process.ROIlabel{1};
                end
            else
                ROIChoice{1} = '';
            end
        end
        
        
        
        if size(Tracks{i},2) < 8
            [ImmTracks, Dummy] = calculateImmobileTracks...
                (Tracks{i}, handles.Parameters.ThreshL,...
                handles.Parameters.minBoundFrames, handles.Parameters.ThreshH, 0);
            
            
            if isempty(ImmTracks)
                TrackLengthHist1{i} = [0 0];
            else
                
                % Identify bound molecules for each of the movies;
                
                % Calculate survival histogram;
                TrackLengthHist1{i} = ...
                    calculateTrackLength(ImmTracks, ...
                    handles.FrameTime1,handles.Parameters.minBoundFrames);
                
                
                
            end
        else
            ImmTracks = Tracks{i}(Tracks{i}(:,2) >= handles.Parameters.minBoundFrames,:);
            ImmTracks = ImmTracks(ImmTracks(:,3) <= handles.Parameters.ThreshL,:);
            if ~isempty(ImmTracks)
                TrackLength = ImmTracks(:,2);
                LongestTrack = max(TrackLength);
                ShortestTrack = handles.Parameters.minBoundFrames;
                
                TrackLengthHist1{i} = zeros(LongestTrack-ShortestTrack + 1,2);
                TrackLengthHist1{i}(:,1) = (ShortestTrack: LongestTrack)* handles.FrameTime1;
                %         TrackLengthHist(:,2) = hist(TrackLength,(ShortestTrack:LongestTrack));
                for  j = 1:length(TrackLength);
                    TrackLengthHist1{i}(1:TrackLength(j)- ShortestTrack + 1,2) = ...
                        TrackLengthHist1{i}(1:TrackLength(j) - ShortestTrack + 1,2)+1;
                end;
            end
        end
        
        
    else
        errordlg('Some of the files have not been preprocessed');
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
for i = 1:nFiles1
    
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
Hist_Matrix = zeros(nMax,nFiles1);

for i = 1:nFiles1
    Hist_Matrix(1:nMax,i) = NParticles{i}(1:nMax,2);
end

CumNParticles1(:,1) = TimePoints * handles.FrameTime1;
CumNParticles1(:,2) = sum(Hist_Matrix, 2);

%Calculate the avg intensity


% Accumulate the different survival histograms;

% find the histogram containing the longest track
nMax = 1;
TimePoints = [];
for i = 1:nFiles1
    
    if length(TrackLengthHist1{i}(:,1)) > nMax;
        
        TimePoints = TrackLengthHist1{i}(:,1);
        nMax = length(TimePoints);
    end
end

if isempty(TimePoints)
    errordlg('Too few bound particles to produce an histogram');
    return
end


% Accumulate the histogram
Hist_Matrix = zeros(nMax,nFiles1);

for i = 1:nFiles1
    lengthHist = length(TrackLengthHist1{i}(:,1));
    Hist_Matrix(1:lengthHist,i) = TrackLengthHist1{i}(:,2);
end

CumTrackLengthHist1(:,1) = TimePoints;
CumTrackLengthHist1(:,2) = sum(Hist_Matrix, 2);

ResTimeHist1(:,1)= CumTrackLengthHist1(:,1);
ResTimeHist1(:,2) = ...
    [(CumTrackLengthHist1(1:end-1,2) - CumTrackLengthHist1(2:end,2));CumTrackLengthHist1(end,2)];
%%%
for i = 1:nFiles2
    
    ROIidx = 0;
    Temp = load([handles.PathName2,handles.FileNames2{i}],'Results');
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
                if size(Temp.Results.Process.ROIlabel,1) > 1
                    if ~isempty(ROIChoice) %if a named ROI was selected before check to see if an ROI with that name exists in this dataset
                        if ~strcmp(ROIChoice{1},'All')
                            for j = 1:length(ROIChoice)
                                for k = 1:length(Temp.Results.Process.ROIlabel)
                                    if strcmpi(ROIChoice{j},Temp.Results.Process.ROIlabel{k})
                                        ROIidx = k;
                                        break
                                    end
                                end
                            end
                        else
                            ROIidx = size(Temp.Results.Process.ROIlabel,1) + 1; % If all was selected before, just set the index to something greater than the number of ROIs
                        end
                    end
                    
                    if ROIidx == 0 %if nothing has been found automatically, provide the user with a dialog box to select the ROI
                        
                        ROIstring = Temp.Results.Process.ROIlabel;
                        ROIstring{end+1,1} = 'All';
                        ROIidx = ROIchooseDlg(ROIstring);
                        
                    end
                    %Update the Tracks & NParticles data
                    if ROIidx <= size(Temp.Results.Process.ROIlabel,1)
                        if size(Tracks{i},2) < 8
                            Tracks{i} = Tracks{i}(Tracks{i}(:,5) == ROIidx,:);
                        else
                            Tracks{i} = Tracks{i}(Tracks{i}(:,9) == ROIidx,:);
                        end
                        NParticles{i} = [NParticles{i}(:,1) NParticles{i}(:,ROIidx+1)];
                        ROIChoice{end+1,1} = Temp.Results.Process.ROIlabel{ROIidx};
                    else
                        tmp = NParticles{i}(:,2:size(Temp.Results.Process.ROIlabel,1)+1);
                        NParticles{i} = [NParticles{i}(:,1) sum(tmp,2)];
                        ROIChoice{end+1,1} = 'All';
                    end
                else
                    ROIChoice{1} = Temp.Results.Process.ROIlabel{1};
                end
            else
                ROIChoice{1} = '';
            end
        end
        
        
        
        if size(Tracks{i},2) < 8
            [ImmTracks, Dummy] = calculateImmobileTracks...
                (Tracks{i}, handles.Parameters.ThreshL,...
                handles.Parameters.minBoundFrames, handles.Parameters.ThreshH, 0);
            
            
            if isempty(ImmTracks)
                TrackLengthHist2{i} = [0 0];
            else
                
                % Identify bound molecules for each of the movies;
                
                % Calculate survival histogram;
                TrackLengthHist2{i} = ...
                    calculateTrackLength(ImmTracks, ...
                    handles.FrameTime2,handles.Parameters.minBoundFrames);
                
                
                
            end
        else
            ImmTracks = Tracks{i}(Tracks{i}(:,2) >= handles.Parameters.minBoundFrames,:);
            ImmTracks = ImmTracks(ImmTracks(:,3) <= handles.Parameters.ThreshL,:);
            if ~isempty(ImmTracks)
                TrackLength = ImmTracks(:,2);
                LongestTrack = max(TrackLength);
                ShortestTrack = handles.Parameters.minBoundFrames;
                
                TrackLengthHist2{i} = zeros(LongestTrack-ShortestTrack + 1,2);
                TrackLengthHist2{i}(:,1) = (ShortestTrack: LongestTrack)* handles.FrameTime2;
                %         TrackLengthHist(:,2) = hist(TrackLength,(ShortestTrack:LongestTrack));
                for  j = 1:length(TrackLength);
                    TrackLengthHist2{i}(1:TrackLength(j)- ShortestTrack + 1,2) = ...
                        TrackLengthHist2{i}(1:TrackLength(j) - ShortestTrack + 1,2)+1;
                end;
            end
        end
        
        
    else
        errordlg('Some of the files have not been preprocessed');
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
for i = 1:nFiles1
    
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
Hist_Matrix = zeros(nMax,nFiles2);

for i = 1:nFiles2
    Hist_Matrix(1:nMax,i) = NParticles{i}(1:nMax,2);
end

CumNParticles2(:,1) = TimePoints * handles.FrameTime2;
CumNParticles2(:,2) = sum(Hist_Matrix, 2);

%Combine the CumNParticles data
NParticles1(:,1) = CumNParticles1(:,1);
NParticles2(:,1) = CumNParticles2(:,1);

NParticles1(:,2) = [CumNParticles1(1:end-1,2) - CumNParticles1(2:end,2) ; CumNParticles1(end,2)];
NParticles2(:,2) = [CumNParticles2(1:end-1,2) - CumNParticles2(2:end,2); CumNParticles2(end,2)];

NParticles = [NParticles1; NParticles2];

NParticles = sortrows(NParticles,1);

dT = diff(NParticles(:,1));

ind = find(dT == 0);
if ~isempty(ind)
    NParticles(ind,2) = NParticles(ind,2) + NParticles(ind+1,2);
end
NParticles(ind+1,:) = [];

CumNParticles(:,1) = NParticles(:,1);
% csum = cumsum(ResTimeHist(:,2));
CumNParticles(:,2) = cumsum(NParticles(:,2),1,'reverse');

%Calculate the avg intensity


% Accumulate the different survival histograms;

% find the histogram containing the longest track
nMax = 1;
TimePoints = [];
for i = 1:nFiles2
    
    if length(TrackLengthHist2{i}(:,1)) > nMax;
        
        TimePoints = TrackLengthHist2{i}(:,1);
        nMax = length(TimePoints);
    end
end

if isempty(TimePoints)
    errordlg('Too few bound particles to produce an histogram');
    return
end


% Accumulate the histogram
Hist_Matrix = zeros(nMax,nFiles1);

for i = 1:nFiles1
    lengthHist = length(TrackLengthHist2{i}(:,1));
    Hist_Matrix(1:lengthHist,i) = TrackLengthHist2{i}(:,2);
end

CumTrackLengthHist2(:,1) = TimePoints;
CumTrackLengthHist2(:,2) = sum(Hist_Matrix, 2);
%%%




% Calculate the fraction of bound molecules.

% First calculate the residence time histogram
ResTimeHist2(:,1)= CumTrackLengthHist2(:,1);
ResTimeHist2(:,2) = ...
    [(CumTrackLengthHist2(1:end-1,2) - CumTrackLengthHist2(2:end,2));CumTrackLengthHist2(end,2)];

ResTimeHist = [ResTimeHist1;ResTimeHist2];
ResTimeHist = sortrows(ResTimeHist,1);

dT = diff(ResTimeHist(:,1));

ind = find(dT == 0);
if ~isempty(ind)
    ResTimeHist(ind,2) = ResTimeHist(ind,2) + ResTimeHist(ind+1,2);
end
ResTimeHist(ind+1,:) = [];

CumTrackLengthHist(:,1) = ResTimeHist(:,1);
% csum = cumsum(ResTimeHist(:,2));
CumTrackLengthHist(:,2) = cumsum(ResTimeHist(:,2),1,'reverse');

% Copy Histogram to handles
handles.CumHist = CumTrackLengthHist;

% Then calculate the total number of bound spots
TotalBoundMolecules = sum(ResTimeHist1(:,2).*ResTimeHist1(:,1)/handles.FrameTime1) + sum(ResTimeHist2(:,2).*ResTimeHist2(:,1)/handles.FrameTime2);

% and the total number of spots


TotalMolecules = sum(CumNParticles(:,2));

% Finally divide the two.
PartialBoundFraction = TotalBoundMolecules/TotalMolecules;

BFerror = (sqrt(TotalBoundMolecules)/TotalBoundMolecules + ...
    sqrt(TotalMolecules)/TotalMolecules)*PartialBoundFraction;


% NOTE: THIS IS ONLY THE PARTIAL BOUND FRACTION
% BECAUSE WE ARE LOOKING ONLY AT PARTICLES BOUND
% FOR MORE THAN Nmin. We HAVE TO FIT THE HISTOGRAM OF RESIDENCE
% TIMES AND EXTRAPOLATE THE TRUE VALUE OF THE BOUND FRACTION

% Photobleaching correction



%If one of the datasets has more time-points than the other, the fitting
%will be weird, so cut off CumNParticles for times > maximum time of the
%shortest dataset

maxT = min(max(CumNParticles1(:,1)),max(CumNParticles2(:,1)));

CumNParticles(CumNParticles(:,1) > maxT,:) = [];

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
    
    binFactor = handles.Parameters.bin + 1;
    nBins = floor(length(ResTimeHist(:,1))/binFactor);
%     if strcmp(answer,'Yes')
        TimeBins = reshape(ResTimeHist(1:nBins*binFactor,1), binFactor, nBins);
        CountsBins = reshape(ResTimeHist(1:nBins*binFactor,2), binFactor, nBins);
%     else
        TimeBins_PB = reshape(ResTimeHist_PB(1:nBins*binFactor,1), binFactor, nBins);
        CountsBins_PB = reshape(ResTimeHist_PB(1:nBins*binFactor,2), binFactor, nBins);
%     end
    
    if binFactor ~= 1
        resTimeHist_Binned = [];
        resTimeHist_Binned(:,1) = mean(TimeBins);
        resTimeHist_Binned(:,2) = sum(CountsBins);
        resTimeHist_Binned_PB = [];
        resTimeHist_Binned_PB(:,1) = mean(TimeBins);
        resTimeHist_Binned_PB(:,2) = sum(CountsBins);
    else
        resTimeHist_Binned = ResTimeHist;
        resTimeHist_Binned_PB = ResTimeHist_PB;
    end
    
    % Remove values < 0 from the histogram
    idx = find(resTimeHist_Binned(:,2) <= 0);
    resTimeHist_Binned2 = resTimeHist_Binned;
    resTimeHist_Binned2(idx,:) = [];
    % Fit the binned residence time histogram with exponentials;
    [fitpar,espSigma, ~]= ExpDecay_fit_resTime(resTimeHist_Binned2, 1 );
    fit = [resTimeHist_Binned(:,1),ExpDecay_fun(fitpar,resTimeHist_Binned(:,1))];
    [fitpar2,espSigma2, ~]= ExpDecay_2Cmp_fit_resTime(resTimeHist_Binned2, [1 0.1]);
    
    Esp_Coef_1 = [fitpar2(1) fitpar2(4)*fitpar2(3)];
    Esp_Fit(:,1) = ExpDecay_fun(Esp_Coef_1,resTimeHist_Binned(:,1));
    
    Esp_Coef_2 = [fitpar2(2) fitpar2(4)*(1 - fitpar2(3))];
    Esp_Fit(:,2) = ExpDecay_fun(Esp_Coef_2,resTimeHist_Binned(:,1));
    
    fit2 = [resTimeHist_Binned(:,1),ExpDecay_2Cmp_fun(fitpar2,resTimeHist_Binned(:,1)),Esp_Fit];
        
    handles.Hist.Res(:, 1:2) = resTimeHist_Binned;
    handles.Hist.Res(:,3) = fit(:,2);
    handles.Hist.Res(:,4) = fit2(:,2);
    handles.Hist.Res(:,5) = fit2(:,3);
    handles.Hist.Res(:,6) = fit2(:,4);
    
    [~, pvalRes, FstatRes] = FtestModelCompare(resTimeHist_Binned(:,2),fit(:,2),fit2(:,2),2,4);

    handles.FitPar.Res = [fitpar, fitpar2, PartialBoundFraction, pvalRes; espSigma, espSigma2, BFerror, FstatRes];
    
    idx = find(resTimeHist_Binned_PB(:,2) <= 0);
    resTimeHist_Binned_PB2 = resTimeHist_Binned_PB;
    resTimeHist_Binned_PB2(idx,:) = [];
    % Fit the binned residence time histogram with exponentials;
    [fitpar_PB,espSigma_PB, ~]= ExpDecay_fit_resTime(resTimeHist_Binned_PB2, 1 );
    fit_PB = [resTimeHist_Binned(:,1),ExpDecay_fun(fitpar_PB,resTimeHist_Binned_PB(:,1))];
    
    [fitpar2_PB,espSigma2_PB, ~]= ExpDecay_2Cmp_fit_resTime(resTimeHist_Binned_PB2, [1 0.1]);
    
    Esp_Coef_1 = [fitpar2_PB(1) fitpar2_PB(4)*fitpar2_PB(3)];
    Esp_Fit(:,1) = ExpDecay_fun(Esp_Coef_1,resTimeHist_Binned_PB(:,1));
    
    Esp_Coef_2 = [fitpar2_PB(2) fitpar2_PB(4)*(1 - fitpar2_PB(3))];
    Esp_Fit(:,2) = ExpDecay_fun(Esp_Coef_2,resTimeHist_Binned_PB(:,1));
    
    fit2_PB = [resTimeHist_Binned_PB(:,1),ExpDecay_2Cmp_fun(fitpar2_PB,resTimeHist_Binned_PB(:,1)),Esp_Fit];
    handles.Hist.Res_PB(:, 1:2) = resTimeHist_Binned_PB;
    handles.Hist.Res_PB(:,3) = fit_PB(:,2);
    handles.Hist.Res_PB(:,4) = fit2_PB(:,2);
    handles.Hist.Res_PB(:,5) = fit2_PB(:,3);
    handles.Hist.Res_PB(:,6) = fit2_PB(:,4);
    [~, pvalRes_PB, FstatRes_PB] = FtestModelCompare(resTimeHist_Binned_PB(:,2),fit_PB(:,2),fit2_PB(:,2),2,4);
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
            ' \pm ', num2str(espSigma(1),2), ' s; C_{eq} = ', num2str(fitpar(2),2)],...
            ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
            ' \pm ', num2str(espSigma2(1),2), 's'], ...
            ['k_{off, 2} = ', num2str(fitpar2(2),3),...
            ' \pm ', num2str(espSigma2(2),2), 's; f_1 = ', num2str(fitpar2(3),3),...
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
    
    
    
    % Fit exponentials to the Survival probability
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
            ' \pm ', num2str(espSigma(1),2), ' s; C_{eq} = ', num2str(fitpar(2),2)],...
            ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
            ' \pm ', num2str(espSigma2(1),2), 's'], ...
            ['k_{off, 2} = ', num2str(fitpar2(2),3),...
            ' \pm ', num2str(espSigma2(2),2), 's; f_1 = ', num2str(fitpar2(3),3),...
            '\pm', num2str(espSigma2(3),2), 'C_{eq} = ',  num2str(fitpar2(4),2)]});
        legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit');
        box on;
        set(handles.FstatText,'String',['F stat: ', num2str(FstatSurv,4)]);
        set(handles.pValText,'String',['p-value: ', num2str(pvalSurv,4)]);
    end
    
    
    
    


figname = get(gcf,'Name');
figname = [figname ' ' ROIChoice{1}];
set(gcf,'Name',figname);



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AnalyzeResTimesDiffRates_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AnalyzeResTimesDiffRates_GUI_OutputFcn(hObject, eventdata, handles)
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
plot(CumNParticles(:,1), CumNParticles(:,2)/max(CumNParticles(:,2)), '+', 'MarkerEdgeColor', [0.5 0.5 0.5]');
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
    DefaultName = handles.PathName1;
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
        ' \pm ', num2str(espSigma(1),2), ' s; C_{eq} = ', num2str(fitpar(2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
        ' \pm ', num2str(espSigma2(1),2), 's'], ...
        ['k_{off, 2} = ', num2str(fitpar2(2),3),...
        ' \pm ', num2str(espSigma2(2),2), 's; f_1 = ', num2str(fitpar2(3),3),...
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
        ' \pm ', num2str(espSigma(1),2), ' s; C_{eq} = ', num2str(fitpar(2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
        ' \pm ', num2str(espSigma2(1),2), 's'], ...
        ['k_{off, 2} = ', num2str(fitpar2(2),3),...
        ' \pm ', num2str(espSigma2(2),2), 's; f_1 = ', num2str(fitpar2(3),3),...
        '\pm', num2str(espSigma2(3),2), 'C_{eq} = ',  num2str(fitpar2(4),2)]});
    legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit');
    box on;
    
    
        
    
    
end
set(handles.FstatText,'String',['F stat: ', num2str(Fstat,4)]);
set(handles.pValText,'String',['p-value: ', num2str(pval,4)]);

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
        ' \pm ', num2str(espSigma(1),2), ' s; C_{eq} = ', num2str(fitpar(2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
        ' \pm ', num2str(espSigma2(1),2), 's'], ...
        ['k_{off, 2} = ', num2str(fitpar2(2),3),...
        ' \pm ', num2str(espSigma2(2),2), 's; f_1 = ', num2str(fitpar2(3),3),...
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
        ' \pm ', num2str(espSigma(1),2), ' s; C_{eq} = ', num2str(fitpar(2),2)],...
        ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
        ' \pm ', num2str(espSigma2(1),2), 's'], ...
        ['k_{off, 2} = ', num2str(fitpar2(2),3),...
        ' \pm ', num2str(espSigma2(2),2), 's; f_1 = ', num2str(fitpar2(3),3),...
        '\pm', num2str(espSigma2(3),2), 'C_{eq} = ',  num2str(fitpar2(4),2)]});
    legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit');
    box on;
    
    
        
    
    
end
set(handles.FstatText,'String',['F stat: ', num2str(Fstat,4)]);
set(handles.pValText,'String',['p-value: ', num2str(pval,4)]);
