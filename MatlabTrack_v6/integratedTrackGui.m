function varargout = integratedTrackGui(varargin)
% INTEGRATEDTRACKGUI M-file for integratedtrackgui.fig
%      INTEGRATEDTRACKGUI, by itself, creates a new INTEGRATEDTRACKGUI or raises the existing
%      singleton*.
%
%      H = INTEGRATEDTRACKGUI returns the handle to a new INTEGRATEDTRACKGUI or the handle to
%      the existing singleton*.
%
%      INTEGRATEDTRACKGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTEGRATEDTRACKGUI.M with the given input arguments.
%
%      INTEGRATEDTRACKGUI('Property','Value',...) creates a new INTEGRATEDTRACKGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before integratedtrackgui_openingfcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to integratedtrackgui_openingfcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help integratedtrackgui

% Last Modified by GUIDE v2.5 24-Apr-2017 13:03:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @TrackGui_OpeningFcn, ...
    'gui_OutputFcn',  @TrackGui_OutputFcn, ...
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


% --- Executes just before integratedtrackgui is made visible.
function TrackGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to integratedtrackgui (see VARARGIN)

% Choose default command line output for integratedtrackgui
handles.output = hObject;

% Get Defaults for Acquisition and analysis parameters.

Defaults = load('IntegratedGUI_defaults.mat', ...
    'ShortestTrack', 'frameTime', 'pixelSize',...
    'JHistBinSize','JHistMaxFrameN', 'MSD_fit_Threshold',...
    'ThreshL','ThreshH', 'minBoundFrames');

% Assign Defaults to handles;
% Parameters regarding the acquisition
handles.Parameters.Acquisition.pixelSize = Defaults.pixelSize;
handles.Parameters.Acquisition.frameTime = Defaults.frameTime;
% Parameters regarding the track checking
handles.Parameters.discardCheck = Defaults.ShortestTrack;

% Parameters regarding the Analysis
handles.Parameters.Analysis.JHistBinSize = Defaults.JHistBinSize;
handles.Parameters.Analysis.JHistMaxFrameN = Defaults.JHistMaxFrameN;
handles.Parameters.Analysis.MSD_fit_Threshold= Defaults.MSD_fit_Threshold;

%Parameters regarding the identification of bound molecules
handles.Parameters.Bound.ThreshL = Defaults.ThreshL;
handles.Parameters.Bound.ThreshH = Defaults.ThreshH;
handles.Parameters.Bound.minBoundFrames = Defaults.minBoundFrames;



% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = TrackGui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on slider movement.
function imageSlider_Callback(hObject, eventdata, handles)
% hObject    handle to imageSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

whatToDraw = get(handles.imagePopup, 'Value');
% Image the selected image
handles.Current.Ix = round(get(hObject,'Value'));
switch whatToDraw
    case 1
        handles.Current.Image = ...
            handles.Data.imageStack(handles.Current.Ix).data;
        %         handles.Current.clims = handles.Data.clims;
    case 2
        handles.Current.Image = ...
            handles.Process.filterStack(handles.Current.Ix).data;
        %         handles.Current.clims = handles.Process.clims;
        %         handles.Current.clims = [0 500];
        
end;
set(handles.axes1,'NextPlot','replacechildren');
imagesc(handles.Current.Image, handles.Current.clims);
% imagesc(handles.Current.Image);
%         axis image;
if get(handles.showParticles, 'Value')
    if handles.isFitPSF
        plotParticle(handles.Tracking.Particles, ...
            handles.Current.Ix, handles.isFitPSF);
    else
        plotParticle(handles.Tracking.Centroids, ...
            handles.Current.Ix, handles.isFitPSF);
    end
end

if get(handles.showTracks, 'Value')
    plotTracks(handles.Current.Tracks, handles.Current.Ix);
end

if isfield(handles.Process,'ROIpos')
    plotROI(handles.Process.ROIpos);
    
end

% Update the image indicator in the panel
set(handles.imageCounter, 'String', ['Image ', num2str(handles.Current.Ix),'/',num2str(handles.Data.nImages)]);

impixelinfo;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function imageSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in imagePopup.
function imagePopup_Callback(hObject, eventdata, handles)
% hObject    handle to imagePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns imagePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imagePopup


whatToDraw = get(hObject, 'Value');
% Image the selected image
switch whatToDraw
    case 1
        handles.Current.Image = ...
            handles.Data.imageStack(handles.Current.Ix).data;
        handles.Current.clims = handles.Data.clims;
        set(handles.BlackValSlider,'Min',0,'Max',handles.Data.MaxBrightnessStack,'SliderStep',[1/(handles.Data.MaxBrightnessStack-1),10/(handles.Data.MaxBrightnessStack-1)]);
        set(handles.WhiteValSlider,'Min',0,'Max',handles.Data.MaxBrightnessStack,'SliderStep',[1/(handles.Data.MaxBrightnessStack-1),10/(handles.Data.MaxBrightnessStack-1)]);
        
        
        
    case 2
        handles.Current.Image = ...
            handles.Process.filterStack(handles.Current.Ix).data;
        handles.Current.clims = handles.Process.clims;
        set(handles.BlackValSlider,'Min',0,'Max',handles.Process.MaxBrightnessStack,'SliderStep',[1/(handles.Process.MaxBrightnessStack-1),10/(handles.Process.MaxBrightnessStack-1)]);
        set(handles.WhiteValSlider,'Min',0,'Max',handles.Process.MaxBrightnessStack,'SliderStep',[1/(handles.Process.MaxBrightnessStack-1),10/(handles.Process.MaxBrightnessStack-1)]);
end
set(handles.axes1,'NextPlot','replacechildren');
imagesc(handles.Current.Image, handles.Current.clims);
% axis image;
set(handles.BlackValEdit,'String',num2str(handles.Current.clims(1)));
set(handles.BlackValSlider,'Value',handles.Current.clims(1));
set(handles.WhiteValEdit,'String',num2str(handles.Current.clims(2)));
set(handles.WhiteValSlider,'Value',handles.Current.clims(2));

if get(handles.showParticles, 'Value')
    
    if handles.isFitPSF
        plotParticle(handles.Tracking.Particles, ...
            handles.Current.Ix, handles.isFitPSF);
    else
        plotParticle(handles.Tracking.Centroids, ...
            handles.Current.Ix, handles.isFitPSF);
    end
end

if isfield(handles.Process,'ROIpos')
    plotROI(handles.Process.ROIpos);
    
end
impixelinfo;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function imagePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imagePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showParticles.
function showParticles_Callback(hObject, eventdata, handles)
% hObject    handle to showParticles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showParticles

isOn = get(hObject, 'Value');

if isOn;
    if handles.isFitPSF
        plotParticle(handles.Tracking.Particles, handles.Current.Ix, ...
            handles.isFitPSF);
    else
        plotParticle(handles.Tracking.Centroids, handles.Current.Ix, ...
            handles.isFitPSF);
    end
else
    delete(findobj(gca,'MarkerSize',12));
end;

% --- Executes on button press in showTracks.
function showTracks_Callback(hObject, eventdata, handles)

isOn = get(hObject, 'Value');
if isOn;
    plotTracks(handles.Current.Tracks, handles.Current.Ix);
else
    delete(findobj(gca,'Color','g'));
end;

function lpEdit_Callback(hObject, eventdata, handles)

function lpEdit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hpEdit_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function hpEdit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in filterButton.
function filterButton_Callback(hObject, eventdata, handles)


%
Handles.Process =[];
Handles.Tracking = [];


% % Disable the options to visualize localized particles
% set(handles.showParticles, 'Value', 0);
% set(handles.showParticles, 'Enable', 'off');
% delete(findobj(gca,'Color','r'));
%
% % Disable the options to visualize the tracks
% set(handles.showTracks, 'Value', 0);
% set(handles.showTracks, 'Enable', 'off');
% delete(findobj(gca,'Color','g'));



nFrames = handles.Data.nImages;
loBP = str2double(get(handles.lpEdit,'String'));
hiBP = str2double(get(handles.hpEdit,'String'));

handles.Parameters.Used.Tracking(1) = loBP;
handles.Parameters.Used.Tracking(2) = hiBP;


for i =  1:nFrames
    set(handles.StatusText,'String',['Filtering Image ', num2str(i), ' of ', num2str(nFrames), '...']);
    drawnow;
    handles.Process.filterStack(i).data = ...
        bpass(handles.Data.imageStack(i).data, loBP, hiBP);
    
    %
    %     bg = [];
    %     for j = 1:size(handles.Process.ROIimage,1)
    %         bg_tmp = handles.Process.filterStack(i).data(handles.Process.ROIimage{j});
    %         bg = [bg; bg_tmp(:)];
    %     end
    %     BG_avg(i,:) = mean(bg);
    %     I_obj(i,:) = calcIntSNR(15,BG_avg(i,:));
    %     handles.Process.filterStack(i).data = imhmax(handles.Process.filterStack(i).data,I_obj(i));
    %
end
set(handles.StatusText,'String','Filtering Images... Done');
drawnow;
handles.Process.clims = [min(min(handles.Process.filterStack(1).data))...
    max(max(handles.Process.filterStack(1).data))];
MaxBrightness = zeros(handles.Data.nImages,1);
for i = 1:handles.Data.nImages
    MaxBrightness(i) = max(handles.Process.filterStack(i).data(:));
end
handles.Process.MaxBrightnessStack = 2*max(MaxBrightness); %give a 100% buffer above the maximum

set(handles.imagePopup, 'Enable', 'on');


if get(handles.imagePopup, 'Value') ==2;
    handles.Current.Image = ...
        handles.Process.filterStack(handles.Current.Ix).data;
    handles.Current.clims = handles.Process.clims;
    set(handles.BlackValSlider,'Min',0,'Max',handles.Process.MaxBrightnessStack,'SliderStep',[1/(handles.Process.MaxBrightnessStack-1),10/(handles.Process.MaxBrightnessStack-1)]);
    set(handles.WhiteValSlider,'Min',0,'Max',handles.Process.MaxBrightnessStack,'SliderStep',[1/(handles.Process.MaxBrightnessStack-1),10/(handles.Process.MaxBrightnessStack-1)]);
    
    set(handles.BlackValEdit,'String',num2str(handles.Current.clims(1)));
    set(handles.BlackValSlider,'Value',handles.Current.clims(1));
    set(handles.WhiteValEdit,'String',num2str(handles.Current.clims(2)));
    set(handles.WhiteValSlider,'Value',handles.Current.clims(2));
    set(handles.axes1,'NextPlot','replacechildren');
    imagesc(handles.Current.Image,handles.Current.clims);
    %     axis image;
end;

%Enable the find Peaks controls
set(handles.findpeakButton,'Enable','on');
set(handles.fitPSF,'Enable','on');
set(handles.thresholdEdit,'Enable','on');
set(handles.windowEdit,'Enable','on');
set(handles.SaveFilterMenu,'Enable','on');

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in loadButton.
function loadButton_Callback(hObject, eventdata, handles)

% Load the image.
[fileName,pathName]  = uigetfile('*.tif','Select stack for read out');
if fileName ~= 0
    set(handles.StatusText,'String','Reading Image...');
    drawnow;
    [imageStack, nImages] = TIFread([pathName, fileName]);
    set(handles.StatusText,'String','Reading Image...Done');
    drawnow;
    % Reset Data fields;
    handles.Data = [];
    
    % Update Handles
    handles.Data.fileName = fileName;
    handles.Data.pathName = pathName;
    handles.Data.imageStack = imageStack;
    handles.Data.nImages = nImages;
    
    % Reset all the fields (opening a new file you will lose the unsaved
    % changes to the current file).
    
    handles.Process = [];
    handles.Process.ROIClass = [];
    handles.Parameters.Used = [];
    handles.Tracking = [];
    handles.Analysis = [];
    handles.PreAnalysis = [];
    handles.Current = [];
    handles.Current.Ix = 1;
    % Select current image and set colormap limits.
    handles.Current.Image = handles.Data.imageStack(handles.Current.Ix).data;
    handles.Data.clims = [min(min(handles.Current.Image))...
        max(max(handles.Current.Image))];
    handles.Current.clims = handles.Data.clims;
    MaxBrightness = zeros(handles.Data.nImages,1);
    for i = 1:handles.Data.nImages
        MaxBrightness(i) = max(handles.Data.imageStack(i).data(:));
    end
    handles.Data.MaxBrightnessStack = 2*max(MaxBrightness); %give a 100% buffer above the maximum
    set(handles.BlackValSlider,'Min',0,'Max',handles.Data.MaxBrightnessStack,'SliderStep',[1/(handles.Data.MaxBrightnessStack-1),10/(handles.Data.MaxBrightnessStack-1)]);
    set(handles.WhiteValSlider,'Min',0,'Max',handles.Data.MaxBrightnessStack,'SliderStep',[1/(handles.Data.MaxBrightnessStack-1),10/(handles.Data.MaxBrightnessStack-1)]);
    
    set(handles.BlackValEdit,'String',num2str(handles.Current.clims(1)));
    set(handles.BlackValSlider,'Value',handles.Current.clims(1));
    set(handles.WhiteValEdit,'String',num2str(handles.Current.clims(2)));
    set(handles.WhiteValSlider,'Value',handles.Current.clims(2));
    
    %Reset the ROI list
    set(handles.ROIList,'Value',1);
    set(handles.ROIList,'String','');
    
    
    % Enable the image slider
    if handles.Data.nImages > 1
        % Set Minimum of image slider
        set(handles.imageSlider,'Min', handles.Current.Ix);
        % Set Maximum of image slider
        set(handles.imageSlider,'Max', handles.Data.nImages);
        % Set Value of image slider
        set(handles.imageSlider,'Value', handles.Current.Ix);
        % Set Step of image slider
        set(handles.imageSlider, 'SliderStep', [1/(handles.Data.nImages-1)...
            1/(handles.Data.nImages-1)]);
        % Set Image counter
        set(handles.imageCounter, 'String', ['Image ', ...
            num2str(handles.Current.Ix),'/',num2str(handles.Data.nImages)]);
    end;
    
    % Disable the options to visualize filtered images
    set(handles.imagePopup, 'Value',1);
    set(handles.imagePopup, 'Enable', 'off');
    
    % Disable the otpions to visualize hand checked tracks
    set(handles.TrackPopUp, 'Value',1);
    set(handles.TrackPopUp, 'Enable', 'off');
    
    % Disable the options to visualize localized particles
    set(handles.showParticles, 'Value', 0);
    set(handles.showParticles, 'Enable', 'off');
    
    % Disable the options to visualize the tracks
    set(handles.showTracks, 'Value', 0);
    set(handles.showTracks, 'Enable', 'off');
    
    %Enable the roi button
    set(handles.roiButton,'Enable','on');
    set(handles.LoadROIfromFile,'Enable','on');
    set(handles.StandardROIbutton,'Enable','on');
    set(handles.RemRefImage,'Enable','off');
    
    %Diable the ROI delete & ROI rename buttons
    set(handles.roiRemovePush,'Enable','off');
    set(handles.roiNamePush,'Enable','off');
    
    %Reset the Classification display
    set(handles.ROIClassCur_text,'String','');
    
    
    %Enable the filter button
    set(handles.filterButton,'Enable','on');
    set(handles.lpEdit,'Enable','on');
    set(handles.hpEdit,'Enable','on');
    
    
    % Display current Image
    
    colormap(gray);
    set(handles.axes1,'NextPlot','replace');
    imagesc(handles.Current.Image, handles.Current.clims);
    axis image;
    impixelinfo;
    set(handles.ImportMMorphROIs,'Enable','on');
end
% Update handles structure
guidata(hObject, handles);



function thresholdEdit_Callback(hObject, eventdata, handles)

function thresholdEdit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function windowEdit_Callback(hObject, eventdata, handles)

function windowEdit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in findpeakbutton.
function findpeakButton_Callback(hObject, eventdata, handles)


% Check if filtered stack is available
if ~isfield(handles.Process,'filterStack')
    msgbox('You need to filter the stack first');
    return
end

if ~isfield(handles.Process,'ROIpos')
    msgbox('You need to add an ROI first');
    return
end


% Verify if tracking data is available;
if ~isempty(handles.Tracking)
    answer = questdlg('Tracking data already available. Do you want to overwrite it?',...
        'Overwrite?','Keep', 'Overwrite', 'Keep');
    switch answer
        case 'Overwrite'
            
            set(handles.showTracks, 'Value', 0);
            set(handles.showTracks, 'Enable', 'off');
            delete(findobj(gca,'Color','g'));
            set(handles.TrackPopUp, 'Value',1);
            set(handles.TrackPopUp, 'Enable', 'off');
            
            % Empty vectors with Tracking and Analysis data;
            handles.Tracking = [];
            handles.PreAnalysis = [];
            handles.Analysis = [];
            
            
            
    end
end









Threshold = str2double(get(handles.thresholdEdit,'String'));
hiBP = str2double(get(handles.hpEdit,'String'));
windowSz = str2double(get(handles.windowEdit,'String'));
handles.isFitPSF = get(handles.fitPSF, 'Value');


% Set Output Parameters
handles.Parameters.Used.Tracking(3) = Threshold;
handles.Parameters.Used.Tracking(4) = windowSz;

set(handles.StatusText,'String','Locating Centroids based on intensity...');
drawnow;
% Calculate Particles centroids
Centroid = findParticles(handles.Process.filterStack, Threshold, hiBP,windowSz);




% Check Particles Inside the ROI and fit PSF
if isfield(handles.Process,'ROIpos')
    set(handles.StatusText,'String','Checking which ROIs Particles are in...');
    drawnow;
    CentroidInRoi = InsideROIcheck2(Centroid, handles.Process.ROIimage);
    
    Centroid = CentroidInRoi;
    
    
    if handles.isFitPSF  % if the 'Fit PSF to peaks' check box is on
        
        % fit PSF to peaks
        set(handles.StatusText,'String','Fitting Centroids to 2D Gaussian...');
        drawnow;
        Particles = peak_fit_psf(handles.Data.imageStack,...
            Centroid,windowSz,windowSz);
        Particles2 = InsideROIcheck2(Particles,handles.Process.ROIimage);
        handles.Tracking.Particles = Particles2;
        plotParticle(Particles2,handles.Current.Ix,handles.isFitPSF);
        
        % Enable filtering based on PSF fit parameters
        %         set(handles.Intensity_thr,'Enable','on');
        %         set(handles.sigmaLow,'Enable','on');
        %         set(handles.sigmaHigh,'Enable','on');
        %         set(handles.filterPeaks,'Enable','on');
        %         set(handles.xPlot,'Enable','on');
        %         set(handles.yPlot,'Enable','on');
        %         set(handles.PlotVariables,'Enable','on');
        
        
        
    else % if the 'Fit PSF to peaks' check box is off
        
        % Plot positions of the centroids
        plotParticle(Centroid,handles.Current.Ix,handles.isFitPSF);
        
        % Disable filtering based on PSF fit parameters
        %         set(handles.Intensity_thr,'Enable','off');
        %         set(handles.sigmaLow,'Enable','off');
        %         set(handles.sigmaHigh,'Enable','off');
        %         set(handles.filterPeaks,'Enable','off');
        %         set(handles.xPlot,'Enable','off');
        %         set(handles.yPlot,'Enable','off');
        %         set(handles.PlotVariables,'Enable','off');
    end
    
else
    if handles.isFitPSF
        set(handles.StatusText,'String','Fitting Centroids to 2D Gaussian...');
        drawnow;
        Particles = peak_fit_psf(handles.Data.imageStack,...
            Centroid,windowSz,windowSz);
        handles.Tracking.Particles = Particles;
        plotParticle(Particles,handles.Current.Ix,handles.isFitPSF);
        
        % Enable filtering based on PSF fit parameters
        %         set(handles.Intensity_thr,'Enable','on');
        %         set(handles.sigmaLow,'Enable','on');
        %         set(handles.sigmaHigh,'Enable','on');
        %         set(handles.filterPeaks,'Enable','on');
        %         set(handles.xPlot,'Enable','on');
        %         set(handles.yPlot,'Enable','on');
        %         set(handles.PlotVariables,'Enable','on');
        
    else
        plotParticle(Centroid,handles.Current.Ix,handles.isFitPSF);
        
        % Disable filtering based on PSF fit parameters
        %         set(handles.Intensity_thr,'Enable','off');
        %         set(handles.sigmaLow,'Enable','off');
        %         set(handles.sigmaHigh,'Enable','off');
        %         set(handles.filterPeaks,'Enable','off');
        %         set(handles.xPlot,'Enable','off');
        %         set(handles.yPlot,'Enable','off');
        %         set(handles.PlotVariables,'Enable','off');
    end
end

set(handles.StatusText,'String','Finding Particles Done');
drawnow;
handles.Tracking.Centroids = Centroid;
handles.Tracking.Peaks = peaks;

set(handles.showParticles, 'Value', 1);
set(handles.showParticles, 'Enable', 'on');

%Enable peak filter &Tracking controls
% set(handles.Intensity_thr,'Enable','on');
% set(handles.sigmaLow,'Enable','on');
% set(handles.sigmaHigh,'Enable','on');
% set(handles.filterPeaks,'Enable','on');
% set(handles.editPeakButton,'Enable','on');
% set(handles.xPlot,'Enable','on');
% set(handles.yPlot,'Enable','on');
% set(handles.PlotVariables,'Enable','on');
set(handles.jumpEdit,'Enable','on');
set(handles.shTrackEdit,'Enable','on');
set(handles.gapsCheck,'Enable','on');
% set(handles.Intensity_thr,'Enable','on');
set(handles.trackButton,'Enable','on');


handles.arePeaksFiltered = 0;


% Update handles structure
guidata(hObject, handles);





% --- Executes during object creation, after setting all properties.
function jumpEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to jumpEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function shTrackEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shTrackEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in trackButton.
function trackButton_Callback(hObject, eventdata, handles)


if isempty (handles.Tracking)
    errordlg('You need to find the position of the particles first!');
    return
end
% Disable the otpions to visualize hand checked tracks

set(handles.TrackPopUp, 'Value',1);
set(handles.TrackPopUp, 'Enable', 'off');

maxJump = str2double(get(handles.jumpEdit,'String'));
Trackparam.mem = str2double(get(handles.gapsCheck, 'String'));
Trackparam.good = str2double(get(handles.shTrackEdit,'String'));
Trackparam.dim         =  2;
Trackparam.quiet       =  0;

handles.Parameters.Used.Tracking(5) = maxJump;
handles.Parameters.Used.Tracking(6) = Trackparam.mem;
handles.Parameters.Used.Tracking(7) = Trackparam.good;

if handles.isFitPSF % if particle position has been evaluated via PSF fitting
    Particles = handles.Tracking.Particles(:,[10 11 6 13]);
else
    Particles = handles.Tracking.Centroids(:,[1 2 6 7]);
end;

% if isfield(handles.Process,'ROIlabel')
%     if size(handles.Process.ROIlabel,1) > 1
%         TrackROIsepAns = questdlg('Do you expect particles to move from one ROI to another?','ROI Particle Tracking','Yes','No','Yes');
%     else
%         TrackROIsepAns = 'Yes';
%     end
% else
TrackROIsepAns = 'Yes';
% end

if strcmp(TrackROIsepAns,'No')
    %testing tracking individual ROIs separately
    
    Tracks = cell(max(Particles(:,4)),1);
    TrkPtsAdded = cell(max(Particles(:,4)),1);
    errorcode = zeros(max(Particles(:,4)),1);
    nTrkPtsAdd = 0;
    ROIstring = handles.Process.ROIlabel;
    for i = 1:max(Particles(:,4))
        Particles2Track = Particles(Particles(:,4) == i,1:3);
        if ~isempty(Particles2Track)
            fprintf('Performing Tracking on %s\n',ROIstring{i,:});
            fprintf('--------------------------\n');
            %             tic;
            [Tracks{i,:}, TrkPtsAdded{i,:}, errorcode(i,:)] = trackfunctIG(Particles2Track,maxJump,Trackparam);
            
            %             [Tracks{i,:}, TrkPtsAdded{i,:}] = ConvexStrict2FrameTracking(Particles2Track);
            %             [Tracks{i,:}, TrkPtsAdded{i,:}] = Convex2FrameTracking(Particles2Track);
            %             [Tracks{i,:}, TrkPtsAdded{i,:}] = ngaTracking(Particles2Track);
            %             tElapsedTracking(i) = toc;
            %             fprintf('Tracking took %2.3d\n',tElapsedTracking(i));
            
            
            nTrkPtsAdd = nTrkPtsAdd + size(TrkPtsAdded{i,:},1);
        else
            fprintf('No particles found in %s, so we cannot track in this ROI\n',ROIstring{i,:});
            Tracks{i,:} = [];
            TrkPtsAdded{i,:} = [];
            errorcode(i,:) = 1;
        end
    end
    %     fprintf('Entire Tracking took %2.3d\n',sum(tElapsedTracking));
    Tracks_all = [];
    Track_ind = 1;
    % TrkPtsAdded2 = [];
    % for i = 1: size(TrkPtsAdded,1)
    %     TrkPtsAdded2 = [TrkPtsAdded2;TrkPtsAdded{i,:}];
    % end
    for i = 1:size(Tracks,1)
        Tracks{i,:}(:,5) = i*ones(size(Tracks{i,:},1),1);
        Tracks_tmp = Tracks{i,:};
        
        for j = 1:max(Tracks_tmp(:,4))
            iTrack = Tracks_tmp(Tracks_tmp(:,4) == j,:);
            if ~isempty(iTrack)
                iTrack(:,4) = Track_ind;
                Track_ind = Track_ind + 1;
                Tracks_all = [Tracks_all;iTrack];
                TrkPtsAdded2{Track_ind,:} = TrkPtsAdded{i,:}{j,:};
            end
        end
        
        %     TrackPtsAdded2{init_TrkPt:fin_TrkPt,:} = TrkPtsAdded{i,:};
        
    end
    
    if ~isempty(Tracks_all)
        Tracks_all2 = sortrows(Tracks_all,3);
        test2 = [];
        used = [];
        t_ind = 1;
        TrkPtsAdded = TrkPtsAdded2;
        TrkPtsAdded2 = cell(size(TrkPtsAdded));
        for i = 1:length(Tracks_all2)
            if isempty(find(Tracks_all2(i,4) == used))
                test = Tracks_all2(Tracks_all2(:,4) == Tracks_all2(i,4),:);
                test(:,4) = t_ind;
                TrkPtsAdded2{t_ind,:} = TrkPtsAdded{Tracks_all2(:,4),:};
                t_ind = t_ind +1;
                test2 = [test2;test];
                
                used = [used; Tracks_all2(i,4)];
            end
        end
        TrkPtsAdded = TrkPtsAdded2;
        Tracks = test2;
    else
        TrkPtsAdded = [];
        Tracks = [];
    end
else
    Particles(Particles(:,1) == 0,:) = [];
    Part_tmp = [];
    for i = 1:max(Particles(:,3))
        PartInCurFrame = Particles(Particles(:,3) == i,:);
        if ~isempty(PartInCurFrame)
            Part_tmp = [Part_tmp; PartInCurFrame];
        else
            addvec = [0 0 i 0];
            Part_tmp = [Part_tmp; addvec];
        end
    end
    Particles = Part_tmp;
    
    set(handles.StatusText,'String','Connecting Particles to Form Tracks...');
    drawnow;
    [Tracks, TrkPtsAdded, errorcode] = trackfunctIG(Particles(:,1:3),maxJump,Trackparam);
    %         [Tracks, TrkPtsAdded] = Convex2FrameTracking(Particles(:,1:3));
    %         [Tracks, TrkPtsAdded] = ngaTracking(Particles(:,1:3));
    %         errorcode = 0;
end




if min(errorcode) == 1;
    errordlg({'Too Few Particles to track.';...
        'Increase the number of particles or allow for shorter tracks'});
elseif min(errorcode) == 2;
    errordlg({'Too Many particles to track.';...
        'Decrease the number of particles or decrease the max jump allowed'});
    
elseif min(errorcode) == 3;
    errordlg({'Someone must have messed with the code';...
        'You are in big trouble'});
    
elseif min(errorcode) == 4;
    errordlg({'It is very difficult for me to track';...
        'especially since the movie you loaded contains only one frame'});
    
elseif min(errorcode) == 5;
    errordlg({'I did not find any particles';...
        'therefore I cannot track'});
else
    if strcmp(TrackROIsepAns,'Yes')
        Tracks = InsideROIcheck2(Tracks,handles.Process.ROIimage); %not
        %     needed if we track in each ROI individually
    end
    
    handles.Tracking.Tracks = Tracks;
    handles.Current.Tracks = handles.Tracking.Tracks;
    
    %Add New Track points to particle list.
    if handles.isFitPSF
        ParticlesNew = handles.Tracking.Particles;
        
        x_ind = 10;
        y_ind = 11;
    else
        ParticlesNew = handles.Tracking.Centroids;
        x_ind = 1;
        y_ind = 2;
    end
    Particles = ParticlesNew;
    
    for i = 1:size(Tracks,1)
        x_pos = Tracks(i,1);
        y_pos = Tracks(i,2);
        frame_num = Tracks(i,3);
        
        pIx1 = find(Particles(:,x_ind) == x_pos & ...
            Particles(:,y_ind) == y_pos & ...
            Particles (:,6) == frame_num);
        if isempty(pIx1)
            ParticleAdd(:,1) = x_pos;
            ParticleAdd(:,2) = y_pos;
            ParticleAdd(:,6) = frame_num;
            if handles.isFitPSF
                ParticleAdd(:,10:11) = ParticleAdd(:,1:2);
                ParticleAdd(:,12) = 0;
                if isfield(handles.Process,'ROIpos')
                    ParticleAdd(:,13) = 0;
                end
            else
                if isfield(handles.Process,'ROIpos')
                    ParticleAdd(:,7) = 0;
                end
            end
            ParticlesNew = [ParticlesNew; ParticleAdd];
        end
        
        
    end
    %     for i = 1: length(TrkPtsAdded)
    %
    %         if i == 147
    %             fh = 10;
    %         end
    %         if ~isempty(TrkPtsAdded{i})
    %             Track_tmp = Tracks(Tracks(:,4) == i,:);
    %             for j = 1:length(TrkPtsAdded{i})
    %                 track_tpoints = Track_tmp(Track_tmp(:,3) == TrkPtsAdded{i}(j),:);
    %                 if ~isempty(track_tpoints)
    %
    %                     %Get Centroid positions before & after gap
    %                     Track_before = Track_tmp(Track_tmp(:,3) == (TrkPtsAdded{i}(j) - 1),:);
    %                     Track_after = Track_tmp(Track_tmp(:,3) == (TrkPtsAdded{i}(j) + 1),:);
    %                     if ~isempty(Track_before) && ~isempty(Track_after)
    %                         pIx1 = find(Particles(:,x_ind) == Track_before(:,1) & ...
    %                             Particles(:,y_ind) == Track_before(:,2) & ...
    %                             Particles (:,6) == Track_before(:,3));
    %
    %                         pIx2 = find(Particles(:,x_ind) == Track_after(:,1) & ...
    %                             Particles(:,y_ind) == Track_after(:,2) & ...
    %                             Particles (:,6) == Track_after(:,3));
    %
    %                         ParticlesAdd(:,1) = (Track_after(:,1) + Track_before(:,1))/2;
    %                         ParticlesAdd(:,2) = (Track_after(:,2) + Track_before(:,2))/2;
    % %                         ParticlesAdd(:,1) = (Particles(pIx2,1) + Particles(pIx1,1))/2;
    % %                         ParticlesAdd(:,2) = (Particles(pIx2,2) + Particles(pIx1,2))/2;
    %                         ParticlesAdd(:,6) = track_tpoints(:,3);
    %                         if handles.isFitPSF
    %                             ParticlesAdd(:,10:11) = track_tpoints(:,1:2);
    %                             ParticlesAdd(:,12) = 0;
    %                             if isfield(handles.Process,'ROIpos')
    %                                 ParticlesAdd(:,13) = 0;
    %                             end
    %                         else
    %                             if isfield(handles.Process,'ROIpos')
    %                                 ParticlesAdd(:,7) = 0;
    %                             end
    %                         end
    %                         ParticlesAdd = InsideROIcheck2(ParticlesAdd,handles.Process.ROIimage,0);
    %                         ParticlesNew = [ParticlesNew; ParticlesAdd];
    %                         if isempty(pIx1)
    %                             ParticleAdd1(:,1) = Track_before(:,1);
    %                             ParticleAdd1(:,2) = Track_before(:,2);
    %                             ParticleAdd1(:,6) = Track_before(:,3);
    %                             if handles.isFitPSF
    %                                 ParticleAdd1(:,10:11) = Track_before(:,1:2);
    %                                 ParticleAdd1(:,12) = 0;
    %                                 if isfield(handles.Process,'ROIpos')
    %                                     ParticleAdd1(:,13) = 0;
    %                                 end
    %                             else
    %                                 if isfield(handles.Process,'ROIpos')
    %                                     ParticleAdd1(:,7) = 0;
    %                                 end
    %                             end
    %                             ParticleAdd1 = InsideROIcheck2(ParticleAdd1,handles.Process.ROIimage,0);
    %                             ParticlesNew = [ParticlesNew; ParticleAdd1];
    %                         end
    %                         if isempty(pIx2)
    %                             ParticleAdd2(:,1) = Track_after(:,1);
    %                             ParticleAdd2(:,2) = Track_after(:,2);
    %                             ParticleAdd2(:,6) = Track_after(:,3);
    %                             if handles.isFitPSF
    %                                 ParticleAdd2(:,10:11) = Track_after(:,1:2);
    %                                 ParticleAdd2(:,12) = 0;
    %                                 if isfield(handles.Process,'ROIpos')
    %                                     ParticleAdd2(:,13) = 0;
    %                                 end
    %                             else
    %                                 if isfield(handles.Process,'ROIpos')
    %                                     ParticleAdd2(:,7) = 0;
    %                                 end
    %                             end
    %                             ParticleAdd2 = InsideROIcheck2(ParticleAdd2,handles.Process.ROIimage,0);
    %                             ParticlesNew = [ParticlesNew; ParticleAdd2];
    %                         end
    %
    %                     end
    %                 end
    %             end
    %         end
    %     end
    ParticlesNew = InsideROIcheck2(ParticlesNew,handles.Process.ROIimage);
    if handles.isFitPSF
        ParticlesNew = sortrows(ParticlesNew,[6,13]);
        handles.Tracking.Particles = ParticlesNew;
    else
        ParticlesNew = sortrows(ParticlesNew,[6,7]);
        handles.Tracking.Centroids = ParticlesNew;
    end
    
    
    
    
    set(handles.showTracks, 'Enable', 'on');
    set(handles.showTracks, 'Value', 1);
    %Enable MAnual checking & analysis
    set(handles.checkButton,'Enable','on');
    set(handles.TrackDimButton,'Enable','on');
    set(handles.Track_Preprocess,'Enable','on');
    set(handles.HistD_Analysis,'Enable','on');
    set(handles.Bound_Analysis,'Enable','on');
    set(handles.Export_Params,'Enable','on');
    
    
    plotTracks(handles.Current.Tracks, handles.Current.Ix);
end
set(handles.StatusText,'String','Connecting Particles to Form Tracks...Done');
drawnow;
% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in roiButton.
function roiButton_Callback(hObject, eventdata, handles)


%Determine state of buttons
showPart_select = get(handles.showParticles,'Value');
showPart_enable = get(handles.showParticles,'Enable');

% Disable the options to visualize localized particles
set(handles.showParticles, 'Value', 0);
set(handles.showParticles, 'Enable', 'off');
% delete(findobj(gca,'Color','r'));

% Disable the ROI buttons
set(handles.roiButton, 'Enable', 'off');
%get the current state of the other ROI buttons
ROIrem_enable = get(handles.roiRemovePush,'Enable');
ROIname_enable = get(handles.roiNamePush,'Enable');

set(handles.roiRemovePush,'Enable','off');
set(handles.roiNamePush,'Enable','off');

TrackPop_enable = get(handles.showTracks,'Enable');
% Disable the otpions to visualize hand checked tracks
% set(handles.TrackPopUp, 'Value',1);
set(handles.TrackPopUp, 'Enable', 'off');

showTrack_select = get(handles.showTracks,'Value');
showTrack_enable = get(handles.showTracks,'Enable');
% Disable the options to visualize the tracks
set(handles.showTracks, 'Value', 0);
set(handles.showTracks, 'Enable', 'off');
delete(findobj(gca,'Color','g'));

% Ask if you want to load a reference image for selecting the ROI
if ~isfield(handles.Process,'RefImage') || isempty(handles.Process.RefImage)
    answer = questdlg('Do you want to open a reference image');
else
    answer = 'Yes';
end

switch answer
    case 'Yes'
        %Open the reference image
        if ~isfield(handles.Process,'RefImage') || isempty(handles.Process.RefImage)
            [fileName, pathName]= ...
                uigetfile ('.tif','Open Reference Image', handles.Data.pathName);
            if fileName == 0
                set(handles.roiButton, 'Enable', 'on');
                return
            else
                [StackRef_stack, nImages] = TIFread([pathName, fileName]);
                RefDesc = 'Sum of Other';
                changeDesc_flag = 1;
            end
            
            
            
            
        else
            StackRef_stack(1).data = handles.Process.RefImage;
            nImages = 1;
            RefDesc = get(handles.RefImageDescription,'String');
            changeDesc_flag = 0;
        end
        
        
    case 'No'
        %Calculate average projection of the images
        if get(handles.imagePopup,'Value') == 1
            StackRef_stack = handles.Data.imageStack;
            nImages = handles.Data.nImages;
            RefDesc = 'Sum of Original';
            changeDesc_flag = 1;
%             SumImage = double(handles.Data.imageStack(1).data);
%             for imageIx =  2:handles.Data.nImages;
%                 SumImage = SumImage +  double(handles.Data.imageStack(imageIx).data);
%             end;
        else
            StackRef_stack = handles.Process.filterStack;
            nImages = handles.Data.nImages;
            RefDesc = 'Sum of Filtered';
            changeDesc_flag = 1;
%             SumImage = double(handles.Process.filterStack(1).data);
%             for imageIx =  2:handles.Data.nImages;
%                 SumImage = SumImage +  double(handles.Process.filterStack(imageIx).data);
%             end;
        end
        
        
    case 'Cancel'
        % Enable the ROI button
        
        set(handles.roiButton, 'Enable', 'on');
        return
        
end
for imgIx = 1:nImages
    StackRef(:,:,imgIx) = double(StackRef_stack(imgIx).data);
end
if ~isfield(handles.Process,'RefImage') || isempty(handles.Process.RefImage)
    answer = questdlg('What type of projection do you want to use for the Reference Image','Projection Type','Sum','Maximum','Sum');
else
    answer = 'Sum';
end

switch answer
    case 'Sum'
        SumImage = sum(StackRef,3);
        
    case 'Maximum'
        SumImage = max(StackRef,[],3);
        if changeDesc_flag == 1;
            RefDesc = ['Max', RefDesc(4:end)];
        end
    
end
   set(handles.RefImageDescription,'String',RefDesc);
        

%Plot average projection of the stack
% if strcmp(answer,'No')
%     if max(max(SumImage))/min(min(SumImage(SumImage > 0))) > 6
%         projlims = [min(min(SumImage(SumImage > 0))) max(max(SumImage))];
%     elseif max(max(SumImage))/min(min(SumImage(SumImage > 0))) > 4
%         projlims = [min(min(SumImage(SumImage > 0))) max(max(SumImage))];
%     else
%         projlims = [min(min(SumImage(SumImage > 0))) max(max(SumImage))];
%     end
% else
projlims = [min(min(SumImage(SumImage > 0))) max(max(SumImage))];
% end

if isfield(handles.Analysis,'HistD')||isfield(handles.Analysis,'BoundF')
    rmDataAns = questdlg('Analysis data already exists that will be removed if you specify a new ROI. Proceed?','Existing Data');
    switch rmDataAns
        case 'Yes'
            handles.Analysis = [];
        case 'No'
            return
    end
    
    
end


%Plot average projection of the stack
set(handles.axes1,'NextPlot','replacechildren');
imagesc(SumImage, projlims);
adjustLims = questdlg('Do you want to adjust the brightness of the reference image','B & C');
while strcmp(adjustLims,'Yes')
    defBW = {num2str(projlims(1)),num2str(projlims(2))};
    num_lines = 1;
    prompt = {'Black','White'};
    black_white = inputdlg(prompt,'Adjust extrema',num_lines,defBW);
    projlims = [str2double(black_white{1}),str2double(black_white{2})];
    set(handles.axes1,'NextPlot','replacechildren');
    imagesc(SumImage, projlims);
    adjustLims = questdlg('Do you want to adjust the brightness of the reference image','B & C');
end
    
% axis image;
if isfield(handles.Process,'ROIpos')
    nROIs = length(handles.Process.ROIpos);
    plotROI(handles.Process.ROIpos);
    
else
    nROIs = 0;
end
set(handles.StatusText,'String','Select the vertices of the ROI');
drawnow;
hPoly = impoly;


ROI_label_def = {['ROI ' num2str(nROIs+1)]};
ROI_label = inputdlg('Create a label for the new ROI:','ROI label',1,ROI_label_def);

if ~isempty(ROI_label)
    ROI_list = get(handles.ROIList,'String');
    if isempty(ROI_list)
        ROI_list = cell(1);
        ROI_list{1,:} = ROI_label{1,:};
    else
        ROI_list{nROIs+1,:} = ROI_label{1,:};
    end
    set(handles.ROIList,'String',ROI_list);
    set(handles.ROIList,'Value',nROIs+1);
    
    handles.Process.ROIpos{nROIs+1,:} = getPosition(hPoly);
    handles.Process.ROIimage{nROIs+1,:} = createMask(hPoly);
    handles.Process.ROIlabel = get(handles.ROIList,'String');
    delete(hPoly);
    set(handles.axes1,'NextPlot','replacechildren');
    imagesc(handles.Current.Image, handles.Current.clims);
    %     axis image;
    plotROI(handles.Process.ROIpos);
    
    handles.Process.RefImage = SumImage;
    set(handles.RemRefImage,'Enable','on');
    %Reassign Particle, Centroid & Tracking data
    if isfield(handles.Tracking,'Particles')
        handles.Tracking.Particles = InsideROIcheck2(handles.Tracking.Particles,handles.Process.ROIimage,0);
    end
    
    if isfield(handles.Tracking,'Centroids')
        handles.Tracking.Centroids = InsideROIcheck2(handles.Tracking.Centroids,handles.Process.ROIimage,0);
    end
    
    if isfield(handles.Tracking,'Tracks')
        handles.Tracking.Tracks = InsideROIcheck2(handles.Tracking.Tracks,handles.Process.ROIimage,0);
    end
    
    if isfield(handles.Tracking,'CheckTracks')
        handles.Tracking.CheckTracks = InsideROIcheck2(handles.Tracking.CheckTracks,handles.Process.ROIimage,0);
    end
    
    if isfield(handles.Tracking,'CheckParticles')
        handles.Tracking.Particles = InsideROIcheck2(handles.Tracking.CheckParticles,handles.Process.ROIimage,0);
    end
    
    if isfield(handles,'PreAnalysis')
        if isfield(handles.PreAnalysis,'Tracks_um')
            pixelSize = handles.Parameters.Acquisition.pixelSize;
            
            fileName = [handles.Data.pathName,handles.Data.fileName];
            if isfield(handles.Tracking,'CheckTracks')
                if size(handles.Tracking.CheckTracks,1) == size(handles.PreAnalysis.Tracks_um,1)
                    handles.PreAnalysis.Tracks_um(:,5) = handles.Tracking.CheckTracks(:,5);
                else
                    [Tracks_um, NParticles, IntensityHist] = preProcess_noGUI(handles.Tracking.CheckTracks,handles.Data.imageStack,...
                        handles.Tracking.CheckParticles, pixelSize, handles.Data.nImages, fileName,handles.Process.ROIpos);
                    handles.PreAnalysis.Tracks_um = Tracks_um;
                    handles.PreAnalysis.NParticles = NParticles;
                    handles.PreAnalysis.IntensityHist = IntensityHist;
                end
            else
                if size(handles.Tracking.Tracks,1) == size(handles.PreAnalysis.Tracks_um,1)
                    handles.PreAnalysis.Tracks_um(:,5) = handles.Tracking.Tracks(:,5);
                else
                    [Tracks_um, NParticles, IntensityHist] = preProcess_noGUI(handles.Tracking.Tracks,handles.Data.imageStack,...
                        handles.Tracking.Particles, pixelSize, handles.Data.nImages, fileName,handles.Process.ROIpos);
                    handles.PreAnalysis.Tracks_um = Tracks_um;
                    handles.PreAnalysis.NParticles = NParticles;
                    handles.PreAnalysis.IntensityHist = IntensityHist;
                end
            end
        end
        if isfield(handles.PreAnalysis,'NParticles')
            nparticles = handles.PreAnalysis.NParticles;
            nparticles(:,end+1) = zeros(size(nparticles,1),1);
            handles.PreAnalysis.NParticles = nparticles;
        end
    end
    % Enable the ROI buttons
    set(handles.roiButton, 'Enable', 'on');
    set(handles.roiRemovePush,'Enable','on');
    set(handles.roiNamePush,'Enable','on');
    set(handles.CopyROI,'Enable','on');
    
    set(handles.showParticles, 'Enable', showPart_enable);
    set(handles.showParticles,'Value',showPart_select);
    
    set(handles.showTracks, 'Enable', showTrack_enable);
    set(handles.showTracks,'Value',showTrack_select);
    
    set(handles.showTracks,'Enable',TrackPop_enable);
    if showPart_select == 1
        if handles.isFitPSF
            plotParticle(handles.Tracking.Particles, handles.Current.Ix, ...
                handles.isFitPSF);
        else
            plotParticle(handles.Tracking.Centroids, handles.Current.Ix, ...
                handles.isFitPSF);
        end
    end
    
    if showTrack_select == 1;
        plotTracks(handles.Current.Tracks, handles.Current.Ix);
    end
    
    %Add classes if desired
    if isempty(handles.Process.ROIClass)
        
        AddClassAns = questdlg('Do you want to separate the ROIs into 2 or more classifications?', 'Classify ROIs','Yes','No','Cancel','No');
        if strcmp(AddClassAns,'Cancel')
            delete(hPoly);
            set(handles.roiButton, 'Enable', 'on');
            set(handles.roiRemovePush,'Enable',ROIrem_enable);
            set(handles.roiNamePush,'Enable',ROIname_enable);
        elseif strcmp(AddClassAns,'No')
            handles.Process.ROIClass = {0};
        else
            ClassName = inputdlg('Specify the name of the new class:','New Classification',1,{'Nucleus'});
            handles.Process.ROIClass = cell(size(handles.Process.ROIpos));
            handles.Process.ROIClass{end,:} = ClassName(1);
            handles.Process.AllROIClasses = ClassName(1);
            set(handles.ROIClassCur_text,'String',ClassName{1});
            set(handles.ChangeROIClass,'Enable','on');
        end
    elseif iscell(handles.Process.ROIClass{1,:})
        ROIstring = handles.Process.AllROIClasses;
        ROIstring{end+1,1} = 'New...';
        ROIidx = ROIClassChooseDlg(ROIstring);
        if ROIidx < size(ROIstring,1)
            handles.Process.ROIClass{end+1,:} = ROIstring{ROIidx,:};
        else
            ClassName = inputdlg('Specify the name of the new class:','New Classification',1,{'Nucleus'});
            %Need to verify that it is really new
            isnew = 1;
            for i = 1:size(handles.Process.AllROIClasses,1)
                if strcmpi(ClassName{1},handles.Process.AllROIClasses{i,:})
                    isnew = 0;
                    break;
                end
            end
            if isnew == 1
                handles.Process.AllROIClasses{end+1,:} = ClassName{1};
                handles.Process.ROIClass{end+1,:} = ClassName{1};
            else
                handles.Process.ROIClass{end+1,:} = handles.Process.AllROIClasses{i,:};
            end
            
            
        end
        set(handles.ROIClassCur_text,'String',handles.Process.ROIClass{end,:});
        set(handles.ChangeROIClass,'Enable','on');
    else
        set(handles.ROIClassCur_text,'String','Not defined');
        set(handles.ChangeROIClass,'Enable','on');
    end
    
else
    delete(hPoly);
    set(handles.roiButton, 'Enable', 'on');
    set(handles.roiRemovePush,'Enable',ROIrem_enable);
    set(handles.roiNamePush,'Enable',ROIname_enable);
    set(handles.axes1,'NextPlot','replacechildren');
    imagesc(handles.Current.Image, handles.Current.clims);
    set(handles.showParticles,'Enable',showPart_enable);
    set(handles.showTracks,'Enable',showTrack_enable);
    if get(handles.showParticles, 'Value')
        if handles.isFitPSF
            plotParticle(handles.Tracking.Particles, ...
                handles.Current.Ix, handles.isFitPSF);
        else
            plotParticle(handles.Tracking.Centroids, ...
                handles.Current.Ix, handles.isFitPSF);
        end
    end
    
    if get(handles.showTracks, 'Value')
        plotTracks(handles.Current.Tracks, handles.Current.Ix);
    end
    
    if isfield(handles.Process,'ROIpos')
        plotROI(handles.Process.ROIpos);
        
    end
    
    
end
set(handles.StatusText,'String','');
drawnow;


% Update handles structure
guidata(hObject, handles);



function Intensity_thr_Callback(hObject, eventdata, handles)
% hObject    handle to Intensity_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Intensity_thr as text
%        str2double(get(hObject,'String')) returns contents of Intensity_thr as a double


% --- Executes during object creation, after setting all properties.
function Intensity_thr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Intensity_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigmaLow_Callback(hObject, eventdata, handles)
% hObject    handle to sigmaLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigmaLow as text
%        str2double(get(hObject,'String')) returns contents of sigmaLow as a double


% --- Executes during object creation, after setting all properties.
function sigmaLow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigmaLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigmaHigh_Callback(hObject, eventdata, handles)
% hObject    handle to sigmaHigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigmaHigh as text
%        str2double(get(hObject,'String')) returns contents of sigmaHigh as a double


% --- Executes during object creation, after setting all properties.
function sigmaHigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigmaHigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in filterPeaks.
function filterPeaks_Callback(hObject, eventdata, handles)

if handles.arePeaksFiltered
    allParticles = handles.Tracking.allParticles;
    
else
    allParticles = handles.Tracking.Particles;
end

filterPkPar(1) = str2double(get(handles.Intensity_thr,'String'));
filterPkPar(2) = str2double(get(handles.sigmaLow,'String'));
filterPkPar(3) = str2double(get(handles.sigmaHigh,'String'));

Particles = filterPeaks(allParticles, filterPkPar);


plotParticle(Particles,handles.Current.Ix, handles.isFitPSF)

handles.Tracking.Particles = Particles;
handles.Tracking.allParticles = allParticles;
handles.arePeaksFiltered = 1;

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in gapsCheck.
function gapsCheck_Callback(hObject, eventdata, handles)
% hObject    handle to gapsCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gapsCheck


% --- Executes on selection change in xPlot.
function xPlot_Callback(hObject, eventdata, handles)
% hObject    handle to xPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns xPlot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xPlot


% --- Executes during object creation, after setting all properties.
function xPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in yPlot.
function yPlot_Callback(hObject, eventdata, handles)
% hObject    handle to yPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns yPlot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from yPlot


% --- Executes during object creation, after setting all properties.
function yPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PlotVariables.
function PlotVariables_Callback(hObject, eventdata, handles)
% hObject    handle to PlotVariables (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


whatX = get(handles.xPlot, 'Value');
whatY = get(handles.yPlot, 'Value');

if handles.arePeaksFiltered
    plotVariables(whatX,whatY, handles.Tracking.allParticles, handles.Tracking.Particles);
else
    plotVariables(whatX,whatY, handles.Tracking.Particles);
end;



function shTrackEdit_Callback(hObject, eventdata, handles)
% hObject    handle to shTrackEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of shTrackEdit as text
%        str2double(get(hObject,'String')) returns contents of shTrackEdit
%        as a double


function jumpEdit_Callback(hObject, eventdata, handles)
% hObject    handle to jumpEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of jumpEdit as text
%        str2double(get(hObject,'String')) returns contents of jumpEdit as
%        a double


% --- Executes on button press in fitPSF.
function fitPSF_Callback(hObject, eventdata, handles)
% hObject    handle to fitPSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fitPSF


% --------------------------------------------------------------------
function SetPar_Callback(hObject, eventdata, handles)
% hObject    handle to SetPar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function SetAcqPar_Callback(hObject, eventdata, handles)
% hObject    handle to SetAcqPar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompt = {'Frame Time [s]','Pixel Size [um]'};
title = 'Set Acquisition Parameters';
def = {num2str(handles.Parameters.Acquisition.frameTime),...
    num2str(handles.Parameters.Acquisition.pixelSize)};

NewPars = inputdlg(prompt,title,1,def);
if ~isempty(NewPars)
    handles.Parameters.Acquisition.pixelSize =str2num(NewPars{2});
    handles.Parameters.Acquisition.frameTime = str2num(NewPars{1});
    
    qstring = {'Do you want to overwrite';...
        'the default acquisition parameters?'};
    Qst = questdlg(qstring,'Overwrite?');
    pixelSize = handles.Parameters.Acquisition.pixelSize;
    if strcmp(Qst, 'Yes')
        frameTime = handles.Parameters.Acquisition.frameTime;
        
        save('IntegratedGUI_defaults.mat', ...
            'frameTime', 'pixelSize','-append');
    end
    
end

% Update handles structure
guidata(hObject, handles);





% --------------------------------------------------------------------
function SetCheckPar_Callback(hObject, eventdata, handles)
% hObject    handle to SetCheckPar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompt = {'Minimum length of tracks to check [frames]'};
title = 'Set Checking Parameters';
def = {num2str(handles.Parameters.discardCheck)};

NewPars = inputdlg(prompt,title,1,def);
if ~isempty(NewPars)
    handles.Parameters.discardCheck = str2num(NewPars{1});
    
    qstring = {'Do you want to overwrite';...
        'the default checking parameters?'};
    Qst = questdlg(qstring,'Overwrite?');
    
    if strcmp(Qst, 'Yes')
        ShortestTrack = handles.Parameters.discardCheck;
        
        save('IntegratedGUI_defaults.mat', ...
            'ShortestTrack','-append');
    end
    
end

% Update handles structure
guidata(hObject, handles);



% --------------------------------------------------------------------
function SetAnPar_Callback(hObject, eventdata, handles)
% hObject    handle to SetAnPar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompt = {'Displacement Histogram Bin Size [um]';...
    'Maximum Frame interval for Displacements Histogram [frames]';...
    'Number of time points for linear fit of MSD'};
title = 'Set Analysis Parameters';
def = {num2str(handles.Parameters.Analysis.JHistBinSize);...
    num2str(handles.Parameters.Analysis.JHistMaxFrameN);...
    num2str(handles.Parameters.Analysis.MSD_fit_Threshold)};


NewPars = inputdlg(prompt,title,1,def);
if ~isempty(NewPars)
    handles.Parameters.Analysis.JHistBinSize = str2num(NewPars{1});
    handles.Parameters.Analysis.JHistMaxFrameN = str2num(NewPars{2});
    handles.Parameters.Analysis.MSD_fit_Threshold = str2num(NewPars{3});
    
    qstring = {'Do you want to overwrite';...
        'the default analysis parameters?'};
    Qst = questdlg(qstring,'Overwrite?');
    
    if strcmp(Qst, 'Yes')
        JHistBinSize = handles.Parameters.Analysis.JHistBinSize;
        JHistMaxFrameN = handles.Parameters.Analysis.JHistMaxFrameN;
        MSD_fit_Threshold = handles.Parameters.Analysis.MSD_fit_Threshold;
        save('IntegratedGUI_defaults.mat', ...
            'JHistBinSize','JHistMaxFrameN','MSD_fit_Threshold','-append');
    end
    
end

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function LoadMat_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Open MatLab file and get variables
[FileName,PathName,FilterIndex] = ...
    uigetfile('.mat','Select Matlab File with tracked data');

if FileName ~=0;
    ListVar = whos('-file',[PathName,FileName]);
    set(handles.StatusText,'String','Loading Tracking data file...');
    drawnow;
else
    return
end
% Reset the data content of the GUI
handles.Process = [];
handles.Parameters.Used = [];
handles.Tracking = [];
handles.Analysis = [];
handles.PreAnalysis = [];
handles.Current = [];
handles.Current.Ix = 1;
handles.isFitPSF = 0;


if ismember('Version',{ListVar.name});
    % THIS IS IF THE MAT FILE HAS BEEN CREATED WITH THE CURRENT VERSION OF
    % THE SOFTWARE.
    IN = load([PathName,FileName],'Version');
    if IN.Version == 2.0;
        
        %Load all the contents of the structure.
        IN = load([PathName,FileName], 'Results');
        
        handles.Data =  IN.Results.Data;
%         handles.Data.fileName = FileName;
%         handles.Data.pathName = PathName;
        Parameters = IN.Results.Parameters;
        handles.Process = IN.Results.Process;
        handles.Tracking = IN.Results.Tracking;
        handles.isFitPSF = IN.Results.isFitPSF;
        handles.PreAnalysis = IN.Results.PreAnalysis;
        if ~isfield(handles.Data,'imageStack')
            try 
                [handles.Data.imageStack,nImages] = TIFread(fullfile(handles.Data.pathName,handles.Data.fileName));
            catch
                [handles.Data.fileName, handles.Data.pathName] = uigetfile('*.tif','Could not open Image, Specify File',handles.Data.pathName);
                if handles.Data.fileName == 0
                    return
                end
                [handles.Data.imageStack,nImages] = TIFread(fullfile(handles.Data.pathName,handles.Data.fileName));
            end
            for i = 1:nImages
                handles.Process.filtStack(i).data = ...
                    bpass(handles.Data.imageStack(i).data,Parameters.Tracking(1), Parameters.Tracking(2));
            end
        end
        
        % Disable some of the controls
        % [They will be re-enabled if the variables are found within the input file]
        
        % Disable the otpions to visualize hand checked tracks
        set(handles.TrackPopUp, 'Value',1);
        set(handles.TrackPopUp, 'Enable', 'off');
        
        % Disable the options to visualize filtered images
        set(handles.imagePopup, 'Value',1);
        set(handles.imagePopup, 'Enable', 'off');
        
        % Disable the options to visualize localized particles
        set(handles.showParticles, 'Value', 0);
        set(handles.showParticles, 'Enable', 'off');
        
        % Disable the options to visualize the tracks
        set(handles.showTracks, 'Value', 0);
        set(handles.showTracks, 'Enable', 'off');
        
        %Enable filter button & roiButton
        set(handles.filterButton,'Enable','on');
        set(handles.roiButton,'Enable','on');
        set(handles.LoadROIfromFile,'Enable','on');
        set(handles.lpEdit,'Enable','on');
        set(handles.hpEdit,'Enable','on');
        if isfield(handles.Process, 'RefImage')
            set(handles.RemRefImage,'Enable','on');
        end
        
        % reset current index.
        handles.Current = [];
        handles.Current.Ix = 1;
        
        % Select current image and set colormap limits.
        handles.Current.Image = handles.Data.imageStack(handles.Current.Ix).data;
        MaxBrightness = zeros(handles.Data.nImages,1);
        MinBrightness = zeros(handles.Data.nImages,1);
        for i = 1:handles.Data.nImages
            MaxBrightness(i) = max(handles.Data.imageStack(i).data(:));
            MinBrightness(i) = min(handles.Data.imageStack(i).data(:));
        end
        handles.Data.MaxBrightnessStack = 2*max(MaxBrightness); %give a 100% buffer above the maximum
        handles.Data.clims = [min(MinBrightness),max(MaxBrightness)];
        handles.Current.clims = handles.Data.clims;
        % Display current Image
        
        colormap(gray);
        set(handles.axes1,'NextPlot','replace');
        imagesc(handles.Current.Image, handles.Current.clims);
        axis image;
        
        % Enable the image slider
        if handles.Data.nImages > 1
            % Set Minimum of image slider
            set(handles.imageSlider,'Min', handles.Current.Ix);
            % Set Maximum of image slider
            set(handles.imageSlider,'Max', handles.Data.nImages);
            % Set Value of image slider
            set(handles.imageSlider,'Value', handles.Current.Ix);
            % Set Step of image slider
            set(handles.imageSlider, 'SliderStep', [1/(handles.Data.nImages-1)...
                1/(handles.Data.nImages-1)]);
            % Set Image counter
            set(handles.imageCounter, 'String', ['Image ', ...
                num2str(handles.Current.Ix),'/',num2str(handles.Data.nImages)]);
        end;
        
        
        % Check if there is the filtered stack
        if isfield(handles.Process, 'filterStack');
            % Enable the options to visualize filtered images
            set(handles.imagePopup, 'Value',1);
            set(handles.imagePopup, 'Enable', 'on');
            set(handles.thresholdEdit,'Enable','on');
            set(handles.windowEdit,'Enable','on');
            set(handles.fitPSF,'Enable','on');
            set(handles.findpeakButton,'Enable','on');
            handles.Process.clims = [min(min(handles.Process.filterStack(handles.Current.Ix).data)) ...
                max(max(handles.Process.filterStack(handles.Current.Ix).data))];
            set(handles.SaveFilterMenu,'Enable','on');
            MaxBrightness = zeros(handles.Data.nImages,1);
            for i = 1:handles.Data.nImages
                MaxBrightness(i) = max(handles.Process.filterStack(i).data(:));
            end
            handles.Process.MaxBrightnessStack = 2*max(MaxBrightness); %give a 100% buffer above the maximum
        end
        %Check if there is ROI data
        if isfield(handles.Process,'ROIpos')
            nROIs = length(handles.Process.ROIpos);
            set(handles.ROIList,'String',handles.Process.ROIlabel);
            plotROI(handles.Process.ROIpos);
            set(handles.ROIList,'Value',1);
            set(handles.roiRemovePush,'Enable','on');
            set(handles.roiNamePush,'Enable','on');
            set(handles.CopyROI,'Enable','on');
            set(handles.StandardROIbutton,'Enable','on');
            set(handles.ChangeROIClass,'Enable','on');
            if ~isfield(handles.Process,'ROIClass')
                handles.Process.ROIClass = [];
            end
            %Reset the Classification display
            
            if ~isempty(handles.Process.ROIClass) && iscell(handles.Process.ROIClass{1,:}(1))
                set(handles.ROIClassCur_text,'String',handles.Process.ROIClass{1,:});
            else
                set(handles.ROIClassCur_text,'String','Not Defined');
            end
        else
            set(handles.ROIClassCur_text,'String','');
            
            
        end
        
        
        
        % Check if there is the data about the particles;
        if isfield(handles.Tracking, 'Centroids');
            % Check if there particles are fit with gaussians
            
            if isfield(handles.Tracking, 'Particles');
                
                Centroids =  handles.Tracking.Particles;
                plotParticle(Centroids,handles.Current.Ix,handles.isFitPSF);
            else
                
                Centroids =  handles.Tracking.Centroids;
                plotParticle(Centroids,handles.Current.Ix,handles.isFitPSF);
            end
            
            set(handles.showParticles, 'Value', 1);
            set(handles.showParticles, 'Enable', 'on');
            handles.arePeaksFiltered = 0;
            %             set(handles.Intensity_thr,'Enable','on');
            %             set(handles.sigmaLow,'Enable','on');
            %             set(handles.sigmaHigh,'Enable','on');
            %             set(handles.filterPeaks,'Enable','on');
            %             set(handles.editPeakButton,'Enable','on');
            %             set(handles.xPlot,'Enable','on');
            %             set(handles.yPlot,'Enable','on');
            %             set(handles.PlotVariables,'Enable','on');
            set(handles.jumpEdit,'Enable','on');
            %             set(handles.sigmaLow,'Enable','on');
            set(handles.shTrackEdit,'Enable','on');
            set(handles.gapsCheck,'Enable','on');
            set(handles.trackButton,'Enable','on');
            set(handles.fitPSF,'Value',handles.isFitPSF);
            
        end
        
        % Check if there is information about tracks;
        if isfield(handles.Tracking, 'Tracks')
            
            handles.Current.Tracks = handles.Tracking.Tracks;
            plotTracks(handles.Current.Tracks, handles.Current.Ix);
            set(handles.showTracks, 'Value', 1);
            set(handles.showTracks, 'Enable', 'on');
            set(handles.checkButton,'Enable','on');
            set(handles.TrackDimButton,'Enable','on');
            set(handles.Track_Preprocess,'Enable','on');
            set(handles.HistD_Analysis,'Enable','on');
            set(handles.Bound_Analysis,'Enable','on');
            set(handles.Export_Params,'Enable','on');
            
        end
        
        
        % Check if there are hand checked tracks
        if isfield(handles.Tracking, 'CheckTracks')
            handles.Current.Tracks = handles.Tracking.CheckTracks;
            set(handles.TrackPopUp, 'Value',2);
            set(handles.TrackPopUp, 'Enable', 'on');
        end
        
        %store tracking parameters;
        if isfield(Parameters, 'Tracking');
            handles.Parameters.Used.Tracking = Parameters.Tracking;
            set(handles.lpEdit, 'String', num2str(handles.Parameters.Used.Tracking(1)));
            set(handles.hpEdit, 'String', num2str(handles.Parameters.Used.Tracking(2)));
            set(handles.thresholdEdit, 'String', num2str(handles.Parameters.Used.Tracking(3)));
            set(handles.windowEdit, 'String', num2str(handles.Parameters.Used.Tracking(4)));
            if length(handles.Parameters.Used.Tracking) > 4
                set(handles.jumpEdit, 'String', num2str(handles.Parameters.Used.Tracking(5)));
                set(handles.gapsCheck, 'String', num2str(handles.Parameters.Used.Tracking(6)));
                set(handles.shTrackEdit, 'String', num2str(handles.Parameters.Used.Tracking(7)));
            end
        end
        
        % check if there are acquisition parameters.
        
        if isfield(Parameters,'Acquisition');
            handles.Parameters.Acquisition = Parameters.Acquisition;
            handles.Parameters.Used.Acquisition = Parameters.Acquisition;
        end
        
        % check if there are checking parameters.
        
        if isfield(Parameters, 'discardCheck');
            handles.Parameters.discardCheck = Parameters.discardCheck;
            handles.Parameters.Used.discardCheck = Parameters.discardCheck;
        end
        
        % check if there are analysis parameters.
        
        if isfield(Parameters, 'Analysis');
            handles.Parameters.Analysis = Parameters.Analysis;
            handles.Parameters.Used.Analysis = Parameters.Analysis;
        end
        
        % check if there are parameters for the bound molecules.
        
        if isfield(Parameters, 'Bound');
            handles.Parameters.Bound = Parameters.Bound;
            handles.Parameters.Used.Bound = Parameters.Bound;
        end
        
        
        
        
        
        
    end
    
    
else
    % THIS IS IF THE MAT FILE HAS BEEN CREATED WITH THE OLD VERSION OF THE
    % SOFTWARE.
    
    % LOAD STACK
    if ismember('Stack', {ListVar.name});
        IN = load([PathName,FileName], 'Stack');
        imageStack = IN.Stack;
        handles.Data.fileName = FileName;
        handles.Data.pathName = PathName;
        handles.Data.imageStack = imageStack;
        handles.Data.nImages = length(imageStack);
        
        
        % Reset all the fields (opening a new file you will lose the unsaved
        % changes to the current file).
        
        handles.Process = [];
        handles.Parameters.Used.Tracking = [];
        handles.Tracking = [];
        handles.Analysis = [];
        handles.PreAnalysis = [];
        handles.Current = [];
        handles.Current.Ix = 1;
        
        % Select current image and set colormap limits.
        handles.Current.Image = handles.Data.imageStack(handles.Current.Ix).data;
        handles.Current.clims = [min(min(handles.Current.Image))...
            max(max(handles.Current.Image))];
        
        
        
        
        
        % Enable the image slider
        if handles.Data.nImages > 1
            % Set Minimum of image slider
            set(handles.imageSlider,'Min', handles.Current.Ix);
            % Set Maximum of image slider
            set(handles.imageSlider,'Max', handles.Data.nImages);
            % Set Value of image slider
            set(handles.imageSlider,'Value', handles.Current.Ix);
            % Set Step of image slider
            set(handles.imageSlider, 'SliderStep', [1/(handles.Data.nImages-1)...
                1/(handles.Data.nImages)]);
            % Set Image counter
            set(handles.imageCounter, 'String', ['Image ', ...
                num2str(handles.Current.Ix),'/',num2str(handles.Data.nImages)]);
        end;
        
        
        % Disable the otpions to visualize hand checked tracks
        set(handles.TrackPopUp, 'Value',1);
        set(handles.TrackPopUp, 'Enable', 'off');
        
        % Disable the options to visualize filtered images
        set(handles.imagePopup, 'Value',1);
        set(handles.imagePopup, 'Enable', 'off');
        
        % Disable the options to visualize localized particles
        set(handles.showParticles, 'Value', 0);
        set(handles.showParticles, 'Enable', 'off');
        
        % Disable the options to visualize the tracks
        set(handles.showTracks, 'Value', 0);
        set(handles.showTracks, 'Enable', 'off');
        
        % Display current Image
        %         axis image;
        colormap(gray);
        set(handles.axes1,'NextPlot','replacechildren');
        imagesc(handles.Current.Image, handles.Current.clims);
        
        
        
    end
    
    % LOAD AND SET TRACKING PARAMETERS
    
    if ismember('AA_Parameters', {ListVar.name});
        
        IN = load([PathName,FileName], 'AA_Parameters');
        handles.Parameters.Used.Tracking = IN.AA_Parameters(1:7);
        
        set(handles.lpEdit, 'String', num2str(IN.AA_Parameters(3)));
        set(handles.hpEdit, 'String', num2str(IN.AA_Parameters(4)));
        set(handles.thresholdEdit, 'String', num2str(IN.AA_Parameters(5)));
        set(handles.windowEdit, 'String', num2str(IN.AA_Parameters(6)));
        set(handles.jumpEdit, 'String', num2str(IN.AA_Parameters(7)));
        set(handles.gapsCheck, 'Value', IN.AA_Parameters(8));
        set(handles.shTrackEdit, 'String', num2str(IN.AA_Parameters(9)));
        
    end
    
    % LOAD CENTROIDS
    if ismember('AA_Centroids', {ListVar.name});
        
        IN = load([PathName,FileName], 'AA_Centroids');
        Centroid = IN.AA_Centroids;
        handles.isFitPSF = 0;
        plotParticle(Centroids,handles.Current.Ix,handles.isFitPSF);
        handles.Tracking.Centroids = Centroid;
        set(handles.showParticles, 'Value', 1);
        set(handles.showParticles, 'Enable', 'on');
        handles.arePeaksFiltered = 0;
    end
    
    
    % LOAD PARTICLES
    if ismember('AA_Particles', {ListVar.name});
        
        IN = load([PathName,FileName], 'AA_Particles');
        Particles = IN.AA_Particles;
        handles.isFitPSF = 1;
        plotParticle(Particles,handles.Current.Ix,handles.isFitPSF);
        handles.Tracking.Particles = Particles;
        set(handles.showParticles, 'Value', 1);
        set(handles.showParticles, 'Enable', 'on');
        handles.arePeaksFiltered = 0;
    end
    
    
    % LOAD TRACKS
    if ismember('BB_pxTracks', {ListVar.name});
        
        IN = load([PathName,FileName], 'BB_pxTracks');
        Tracks = IN.BB_pxTracks;
        handles.Tracking.Tracks = Tracks;
        handles.Current.Tracks = handles.Tracking.Tracks;
        plotTracks(handles.Current.Tracks, handles.Current.Ix);
        set(handles.showTracks, 'Value', 1);
        set(handles.showTracks, 'Enable', 'on');
        handles.arePeaksFiltered = 0;
        handles.Current.Tracks = handles.Tracking.Tracks;
        set(handles.TrackPopUp, 'Value', 1);
        set(handles.TrackPopUp, 'Enable', 'off');
    end
    
    % LOAD EDITED TRACKS
    
    if ismember('BB_pxTr_edit', {ListVar.name});
        
        IN = load([PathName,FileName], 'BB_pxTr_edit');
        Tracks = IN.BB_pxTr_edit;
        handles.Tracking.CheckTracks = Tracks;
        set(handles.showTracks, 'Value', 1);
        set(handles.showTracks, 'Enable', 'on');
        
        
        % Enable the otpions to visualize hand checked tracks
        set(handles.TrackPopUp, 'Value', 2);
        set(handles.TrackPopUp, 'Enable', 'on');
        handles.Current.Tracks = handles.Tracking.CheckTracks;
        plotTracks(handles.Current.Tracks,handles.Current.Ix);
        
    end
    
end
set(handles.StatusText,'String','Loading Tracking data file...Done');
drawnow;
% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function SaveMat_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Assign .mat file where variables will be saved.
defName = [handles.Data.pathName, handles.Data.fileName];
SearchStr = '(.*)\.\w*';
defName = regexprep(defName, SearchStr, '$1');
FilterSpec = {'*.mat'};
[FileNameOut,PathNameOut] = uiputfile(FilterSpec,'Save Variables',defName);
if FileNameOut ~=0
    set(handles.StatusText,'String','Saving MAT file...');
    drawnow;
    Version = 2.0;
    
    %     save([PathNameOut, FileNameOut], 'Version');
    Results.Data = handles.Data;
    Results.Parameters = handles.Parameters.Used;
    Results.Process = handles.Process;
    Results.Tracking = handles.Tracking;
    Results.Analysis = handles.Analysis;
    Results.PreAnalysis = handles.PreAnalysis;
    Results.isFitPSF = handles.isFitPSF;
    save([PathNameOut, FileNameOut], 'Version','Results', '-v7.3');
    set(handles.StatusText,'String','Saving MAT file...Done');
    drawnow;
end


% --------------------------------------------------------------------
function ExportToWorkSpace_Callback(hObject, eventdata, handles)
% hObject    handle to ExportToWorkSpace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

assignin('base','AA_Centroids',handles.Tracking.Centroids);
if handles.isFitPSF % if particle position has been evaluated via PSF fitting
    assignin('base','AA_Particles',handles.Tracking.Particles);
end

assignin('base','BB_pxTracks',handles.Tracking.Tracks);
assignin('base','BB_Tracking_Parameters', handles.Parameters.Used.Tracking);
assignin('base','Stack', handles.Data.imageStack);
if isfield(handles.Tracking,'CheckTracks')
    assignin('base','BB_pxTr_Edit', handles.Tracking.CheckTracks);
end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function checkButton_Callback(hObject, eventdata, handles)

handles.PreAnalysis = [];
handles.Analysis = [];
handles.Parameters.Used.discardCheck = handles.Parameters.discardCheck;

if ~isfield(handles.Tracking, 'Tracks')
    errordlg('You need to track the particles first!');
    return
end;



if isfield(handles.Tracking, 'CheckTracks')
    
    button = questdlg(['Some tracks have already been edited. ';...
        'Do you want to continue or start over?'],'Edited Tracks Found', ...
        'Continue', 'Start Over','Continue');
    
else
    button = 'Start Over';
    
end

switch button
    
    case 'Start Over'
        TracksTemp = discardTracks(handles.Tracking.Tracks, ...
            handles.Parameters.discardCheck);
        if ~isfield(handles,'Process')
            handles.Process.ROIpos = [];
        else
            if ~isfield(handles.Process,'ROIpos')
                handles.Process.ROIpos = [];
            end
        end
        %         AnalysisType = questdlg('Do you want to adjust the individual Track Points, or just define the Track length & width?','Type of Checking','Points','Dimensions','Points');
        AnalysisType = 'Points';
        if handles.isFitPSF
            if strcmp(AnalysisType,'Points')
                hCheckGUI = checkTracks2_GUI...
                    (handles.Data.imageStack, TracksTemp, handles.Tracking.Particles, handles.Process);
                %                 [handles.Tracking.CheckTracks handles.Tracking.CheckParticles] = checkTracks2_GUI...
                %                     (handles.Data.imageStack, TracksTemp, handles.Tracking.Particles, handles.Process);
                waitfor(hCheckGUI);
                handles.Tracking.CheckTracks = getappdata(0,'outTracks');
                handles.Tracking.CheckParticles = getappdata(0,'outParticles');
                handles.Tracking.TrackDim = [];
            else
                handles.Tracking.TrackDim = zeros(max(TracksTemp(:,4)),9);
                [handles.Tracking.CheckTracks handles.Tracking.CheckParticles,handles.Tracking.TrackDim] = checkTracks3_GUI...
                    (handles.Data.imageStack, TracksTemp, handles.Tracking.Particles, handles.Process,handles.Tracking.TrackDim);
                handles.Tracking.TrackDim(handles.Tracking.TrackDim(:,1) == 0,:) = [];
            end
        else
            if strcmp(AnalysisType,'Points')
                hCheckGUI = checkTracks2_GUI...
                    (handles.Data.imageStack, TracksTemp, handles.Tracking.Centroids, handles.Process);
                %                 [handles.Tracking.CheckTracks handles.Tracking.CheckParticles] = checkTracks2_GUI...
                %                     (handles.Data.imageStack, TracksTemp, handles.Tracking.Centroids, handles.Process);
                waitfor(hCheckGUI);
                handles.Tracking.CheckTracks = getappdata(0,'outTracks');
                handles.Tracking.CheckParticles = getappdata(0,'outParticles');
                handles.Tracking.TrackDim = [];
            else
                handles.Tracking.TrackDim = zeros(max(TracksTemp(:,4)),9);
                [handles.Tracking.CheckTracks handles.Tracking.CheckParticles,handles.Tracking.TrackDim] = checkTracks3_GUI...
                    (handles.Data.imageStack, TracksTemp, handles.Tracking.Centroids, handles.Process,handles.Tracking.TrackDim);
                handles.Tracking.TrackDim(handles.Tracking.TrackDim(:,1) == 0,:) = [];
            end
        end
        
    case 'Continue'
        %         if ~isempty(handles.Tracking.TrackDim)
        %             AnalysisType = questdlg('Do you want to adjust the individual Track Points, or just define the Track length & width?','Type of Checking','Points','Dimensions','Points');
        %         else
        AnalysisType = 'Points';
        %         end
        if isfield(handles.Tracking, 'CheckParticles')
            if strcmp(AnalysisType,'Points')
                hCheckGUI = checkTracks2_GUI...
                    (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.CheckParticles, handles.Process);
                %                 [handles.Tracking.CheckTracks handles.Tracking.CheckParticles] = checkTracks2_GUI...
                %                     (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.CheckParticles, handles.Process);
                waitfor(hCheckGUI);
                handles.Tracking.CheckTracks = getappdata(0,'outTracks');
                handles.Tracking.CheckParticles = getappdata(0,'outParticles');
                handles.Tracking.TrackDim = [];
            else
                [handles.Tracking.CheckTracks handles.Tracking.CheckParticles,handles.Tracking.TrackDim] = checkTracks3_GUI...
                    (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.CheckParticles, handles.Process,handles.Tracking.TrackDim);
                handles.Tracking.TrackDim(handles.Tracking.TrackDim(:,1) == 0,:) = [];
            end
        else
            
            if handles.isFitPSF
                if strcmp(AnalysisType,'Points')
                    hCheckGUI = checkTracks2_GUI...
                        (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.Particles, handles.Process);
                    %                     [handles.Tracking.CheckTracks handles.Tracking.CheckParticles] = checkTracks2_GUI...
                    %                         (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.Particles, handles.Process);
                    waitfor(hCheckGUI);
                    handles.Tracking.CheckTracks = getappdata(0,'outTracks');
                    handles.Tracking.CheckParticles = getappdata(0,'outParticles');
                    handles.Tracking.TrackDim = [];
                else
                    [handles.Tracking.CheckTracks handles.Tracking.CheckParticles,handles.Tracking.TrackDim] = checkTracks3_GUI...
                        (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.Particles, handles.Process,handles.Tracking.TrackDim);
                    handles.Tracking.TrackDim(handles.Tracking.TrackDim(:,1) == 0,:) = [];
                    
                end
            else
                if strcmp(AnalysisType,'Points')
                    hCheckGUI = checkTracks2_GUI...
                        (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.Centroids, handles.Process);
                    %                     [handles.Tracking.CheckTracks handles.Tracking.CheckParticles] = checkTracks2_GUI...
                    %                         (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.Centroids, handles.Process);
                    waitfor(hCheckGUI);
                    handles.Tracking.CheckTracks = getappdata(0,'outTracks');
                    handles.Tracking.CheckParticles = getappdata(0,'outParticles');
                    
                    handles.Tracking.TrackDim = [];
                else
                    [handles.Tracking.CheckTracks handles.Tracking.CheckParticles,handles.Tracking.TrackDim] = checkTracks3_GUI...
                        (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.Centroids, handles.Process,handles.Tracking.TrackDim);
                    handles.Tracking.TrackDim(handles.Tracking.TrackDim(:,1) == 0,:) = [];
                end
            end
        end
end



% Enable the otpions to visualize hand checked tracks
set(handles.TrackPopUp, 'Value', 2);
set(handles.TrackPopUp, 'Enable', 'on');
handles.Current.Tracks = handles.Tracking.CheckTracks;
plotTracks(handles.Current.Tracks,handles.Current.Ix);


% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in TrackPopUp.
function TrackPopUp_Callback(hObject, eventdata, handles)
% hObject    handle to TrackPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns TrackPopUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TrackPopUp

whatToDraw = get(hObject, 'Value');
% Image the selected image
switch whatToDraw
    case 1
        handles.Current.Tracks = handles.Tracking.Tracks;
        plotTracks(handles.Current.Tracks,handles.Current.Ix);
    case 2
        handles.Current.Tracks = handles.Tracking.CheckTracks;
        plotTracks(handles.Current.Tracks,handles.Current.Ix);
        
        %         axis image;
end;

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function TrackPopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrackPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function SetBoundPar_Callback(hObject, eventdata, handles)
% hObject    handle to SetBoundPar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


prompt = {'Frame to Frame max displacement for bound molecules [um]';...
    'Maximum displacement for bound molecules[um]';...
    'Minimum number of frames to call the molecule bound [frames]'};
title = 'Set Parameters for the isolation of bound molecules';
def = {num2str(handles.Parameters.Bound.ThreshL);...
    num2str(handles.Parameters.Bound.ThreshH);...
    num2str(handles.Parameters.Bound.minBoundFrames)};


NewPars = inputdlg(prompt,title,1,def);
if ~isempty(NewPars)
    handles.Parameters.Bound.ThreshL = str2num(NewPars{1});
    handles.Parameters.Bound.ThreshH = str2num(NewPars{2});
    handles.Parameters.Bound.minBoundFrames = str2num(NewPars{3});
    
    qstring = {'Do you want to overwrite';...
        'the default parameters for the selection of bound molecules?'};
    Qst = questdlg(qstring,'Overwrite?');
    
    if strcmp(Qst, 'Yes')
        ThreshL = handles.Parameters.Bound.ThreshL;
        ThreshH = handles.Parameters.Bound.ThreshH;
        minBoundFrames = handles.Parameters.Bound.minBoundFrames;
        save('IntegratedGUI_defaults.mat', ...
            'ThreshL','ThreshH','minBoundFrames','-append');
    end
    
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in Track_Preprocess.
function Track_Preprocess_Callback(hObject, eventdata, handles)
% hObject    handle to Track_Preprocess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get Tracks and send them to the preprocess GUI
handles.Analysis = [];
handles.Parameters.Used.Acquisition = handles.Parameters.Acquisition;
if ~isfield(handles.Tracking,'Tracks');
    errordlg('You need to track the molecules first!','Error')
    
else
    set(handles.StatusText,'String','Pre-processing Tracks...');
    drawnow;
    if isfield(handles.Tracking, 'CheckTracks')
        answer = questdlg('Do you want to analyze the hand-checked tracks');
        switch answer
            case 'Yes'
                if ~isempty(handles.Tracking.TrackDim)
                    whichTracks = questdlg('Do you want to process the track points, or dimensions?','What to Process','Points','Dimensions','Dimensions');
                else
                    whichTracks = 'Points';
                end
                if strcmp(whichTracks,'Points');
                    Tracks = handles.Tracking.CheckTracks;
                else
                    Tracks = handles.Tracking.TrackDim;
                end
                if isfield (handles.Tracking, 'CheckParticles')
                    Particles = handles.Tracking.CheckParticles;
                else
                    if handles.isFitPSF
                        Particles = handles.Tracking.Particles;
                    else
                        Particles = handles.Tracking.Centroids;
                    end
                end
            case 'No'
                Tracks = handles.Tracking.Tracks;
                if handles.isFitPSF
                    Particles = handles.Tracking.Particles;
                else
                    Particles = handles.Tracking.Centroids;
                end
        end
    else
        Tracks = handles.Tracking.Tracks;
        if handles.isFitPSF
            Particles = handles.Tracking.Particles;
        else
            Particles = handles.Tracking.Centroids;
        end
    end
    
    pixelSize = handles.Parameters.Acquisition.pixelSize;
    
    fileName = [handles.Data.pathName,handles.Data.fileName];
    
    
    %     if handles.isFitPSF
    [Tracks_um, NParticles, IntensityHist] = preProcess_GUI(Tracks,handles.Data.imageStack,...
        Particles, pixelSize, handles.Data.nImages, fileName,handles.Process.ROIpos);
    %     else
    %         [Tracks_um, NParticles, IntensityHist] = preProcess_GUI(Tracks,handles.Data.imageStack,...
    %             handles.Tracking.Centroids, pixelSize, handles.Data.nImages, fileName,handles.Process.ROIpos);
    %     end
    
    handles.PreAnalysis.Tracks_um = Tracks_um;
    handles.PreAnalysis.NParticles = NParticles;
    handles.PreAnalysis.IntensityHist = IntensityHist;
    set(handles.StatusText,'String','Pre-processing Tracks...Done');
    drawnow;
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in HistD_Analysis.
function HistD_Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to HistD_Analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.PreAnalysis)
    errordlg('You need to preprocess the tracks first!','Error!');
else
    Tracks_all = handles.PreAnalysis.Tracks_um;
    if isfield(handles.Process,'ROIpos')
        if length(handles.Process.ROIpos) > 1
            ROIanswer = questdlg('Analyze displacements in selected ROI or all?', 'ROI analysis','ROI','All','Cancel','ROI');
            switch ROIanswer
                case 'ROI'
                    ROIval = get(handles.ROIList,'Value');
                    ROInames = get(handles.ROIList,'String');
                    ROIname = ROInames{ROIval,:};
                    
                    Tracks = Tracks_all(Tracks_all(:,5) == ROIval,:);
                case 'All'
                    ROIname = 'All';
                    Tracks = Tracks_all;
                    ROIval = 0;
                case 'Cancel'
                    return
            end
        else
            Tracks = Tracks_all;
            ROIname = '';
            ROIval = 1;
        end
    else
        Tracks = Tracks_all;
        ROIname = '';
        ROIval = -1;
    end
    handles.Parameters.Used.Analysis = handles.Parameters.Analysis;
    fileName = [handles.Data.pathName,handles.Data.fileName];
    
    Parameters = handles.Parameters.Analysis;
    FrameTime =  handles.Parameters.Acquisition.frameTime;
    MaxJump = handles.Parameters.Acquisition.pixelSize *...
        handles.Parameters.Used.Tracking(5);
    
    [HistD, HistD1Fit, HistD1Coefs] = ...
        JDH1_Analysis_GUI(Tracks, Parameters, FrameTime, MaxJump, fileName, ROIname);
    
    if ROIval == -1
        
        handles.Analysis.HistD = HistD;
        handles.Analysis.HistD1Fit = HistD1Fit;
        handles.Analysis.HistD1Coefs = HistD1Coefs;
    elseif ROIval == 0
        nROIs = length(handles.Process.ROIpos);
        if ~isfield(handles.Analysis,'HistD')
            handles.Analysis.HistD = cell(nROIs+1,1);
            handles.Analysis.HistD1Fit = cell(nROIs+1,1);
            handles.Analysis.HistD1Coefs = cell(nROIs+1,1);
        end
        handles.Analysis.HistD{nROIs+1,:} = HistD;
        handles.Analysis.HistD1Fit{nROIs+1,:} = HistD1Fit;
        handles.Analysis.HistD1Coefs{nROIs+1,:} = HistD1Coefs;
    else
        nROIs = length(handles.Process.ROIpos);
        if ~isfield(handles.Analysis,'HistD')
            handles.Analysis.HistD = cell(nROIs+1,1);
            handles.Analysis.HistD1Fit = cell(nROIs+1,1);
            handles.Analysis.HistD1Coefs = cell(nROIs+1,1);
        end
        handles.Analysis.HistD{ROIval,:} = HistD;
        handles.Analysis.HistD1Fit{ROIval,:} = HistD1Fit;
        handles.Analysis.HistD1Coefs{ROIval,:} = HistD1Coefs;
    end
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in MSD_Analysis.
function MSD_Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to MSD_Analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Bound_Analysis.
function Bound_Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to Bound_Analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.PreAnalysis)
    errordlg('You need to preprocess the tracks first!','Error!');
    return
end
Tracks_all = handles.PreAnalysis.Tracks_um;
if isfield(handles.Process,'ROIpos')
    if length(handles.Process.ROIpos) > 1
        ROIanswer = questdlg('Analyze displacements in selected ROI or all?', 'ROI analysis','ROI','All','Cancel','ROI');
        switch ROIanswer
            case 'ROI'
                ROIval = get(handles.ROIList,'Value');
                ROInames = get(handles.ROIList,'String');
                ROIname = ROInames{ROIval,:};
                
                Tracks = Tracks_all(Tracks_all(:,5) == ROIval,:);
                NParticles = [handles.PreAnalysis.NParticles(:,1) handles.PreAnalysis.NParticles(:,ROIval+1)];
            case 'All'
                ROIname = 'All';
                Tracks = Tracks_all;
                ROIval = 0;
                NParticles = [handles.PreAnalysis.NParticles(:,1) sum(handles.PreAnalysis.NParticles(:,2:end),2)];
            case 'Cancel'
                return
        end
    else
        Tracks = Tracks_all;
        ROIname = '';
        ROIval = 1;
        NParticles = [handles.PreAnalysis.NParticles(:,1) sum(handles.PreAnalysis.NParticles(:,2:end),2)];
    end
else
    Tracks = Tracks_all;
    ROIname = '';
    ROIval = -1;
    NParticles = [handles.PreAnalysis.NParticles(:,1) sum(handles.PreAnalysis.NParticles(:,2:end),2)];
end
fileName = [handles.Data.pathName,handles.Data.fileName];
handles.Parameters.Used.Bound = handles.Parameters.Bound;

Parameters = handles.Parameters.Bound;
FrameTime =  handles.Parameters.Acquisition.frameTime;

nImages = handles.Data.nImages;


[BoundF, ResTimeH,BoundTracks] = ...
    BoundMol_GUI(Tracks, Parameters, FrameTime, NParticles,nImages,fileName,ROIname);



if ROIval == -1
    
    handles.Analysis.BoundF = BoundF;
    handles.Analysis.ResTimeH = ResTimeH;
    handles.Analysis.BoundTracks = BoundTracks;
elseif ROIval == 0
    nROIs = length(handles.Process.ROIpos);
    if ~isfield(handles.Analysis,'BoundF')
        handles.Analysis.BoundF = cell(nROIs+1,1);
        handles.Analysis.ResTimeH = cell(nROIs+1,1);
        handles.Analysis.BoundTracks = cell(nROIs+1,1);
    end
    handles.Analysis.BoundF{nROIs+1,:} = BoundF;
    handles.Analysis.ResTimeH{nROIs+1,:} = ResTimeH;
    handles.Analysis.BoundTracks{nROIs+1,:} = BoundTracks;
    
else
    nROIs = length(handles.Process.ROIpos);
    if ~isfield(handles.Analysis,'BoundF')
        handles.Analysis.BoundF = cell(nROIs+1,1);
        handles.Analysis.ResTimeH = cell(nROIs+1,1);
        handles.Analysis.BoundTracks = cell(nROIs+1,1);
    end
    handles.Analysis.BoundF{ROIval,:} = BoundF;
    handles.Analysis.ResTimeH{ROIval,:} = ResTimeH;
    handles.Analysis.BoundTracks{ROIval,:} = BoundTracks;
    
end



% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in Export_Params.
function Export_Params_Callback(hObject, eventdata, handles)


% Create a Vector with the analysis Params and export it to clipboard

% Copy Tracking Parameters
if ~isfield (handles.Parameters.Used, 'Tracking');
    errordlg('You need to track your data first!');
    return
end
Params = handles.Parameters.Used.Tracking;

if isfield (handles,'isFitPSF')
    Params(8) = handles.isFitPSF;
end

% Copy Check parameters
if isfield (handles.Parameters.Used, 'discardCheck')
    Params(9) = handles.Parameters.Used.discardCheck;
end

% Copy Acquisition parameters
if isfield (handles.Parameters.Used,'Acquisition')
    Params(10:11) = [handles.Parameters.Used.Acquisition.pixelSize,...
        handles.Parameters.Used.Acquisition.frameTime];
end

% Copy Analysis Parameters
if isfield (handles.Parameters.Used,'Analysis')
    Params(12:14) = [handles.Parameters.Used.Analysis.JHistBinSize,...
        handles.Parameters.Used.Analysis.JHistMaxFrameN,...
        handles.Parameters.Used.Analysis.MSD_fit_Threshold];
end

% Copy Analysis Parameters
if isfield (handles.Parameters.Used,'Bound')
    Params(15:17) = [handles.Parameters.Used.Bound.ThreshL...
        handles.Parameters.Used.Bound.ThreshH,...
        handles.Parameters.Used.Bound.minBoundFrames];
end

num2clip(Params);
msgbox('The parameters have been copied to the clipboard');










% hObject    handle to Export_Params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MoreAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to MoreAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Tracks_Col_Callback(hObject, eventdata, handles)
% hObject    handle to Tracks_Col (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if ~isfield(handles.Tracking,'Tracks');
    errordlg('You need to track the molecules first!','Error')
    return;
end;

if isfield(handles.Tracking, 'CheckTracks')
    answer = questdlg('Do you want to analyze the hand-checked tracks');
    switch answer
        case 'Yes'
            Tracks = handles.Tracking.CheckTracks;
        case 'No'
            Tracks = handles.Tracking.Tracks;
    end
end

frameTime = handles.Parameters.Acquisition.frameTime;
pixelSize = handles.Parameters.Acquisition.pixelSize;
PathName = handles.Data.pathName;

if isstruct(handles.Analysis) && isfield(handles.Analysis,'BoundTracks');
    answer = questdlg('Do you want to analyze the track segments identified as bound?');
    switch answer
        case 'Yes'
            pixelSize = handles.Parameters.Used.Acquisition.pixelSize;
            Tracks = handles.Analysis.BoundTracks;
            Tracks(:,1:2)=Tracks(:,1:2)/pixelSize;
    end
end




frameTime = handles.Parameters.Acquisition.frameTime;
pixelSize = handles.Parameters.Acquisition.pixelSize;
PathName = handles.Data.pathName;

ColocalizeTracksGUI(Tracks, pixelSize, frameTime, PathName);


% --------------------------------------------------------------------
function ResTime_Merge_Callback(hObject, eventdata, handles)
% hObject    handle to ResTime_Merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Select multiple Files to open;

if isfield(handles, 'Data')
    [FileNames,PathName,FilterIndex] = uigetfile('.mat',...
        'Select Mat files for the analysis of the residence time',...
        [handles.Data.pathName, handles.Data.fileName],'MultiSelect','on');
else
    [FileNames,PathName,FilterIndex] = uigetfile('.mat',...
        'Select Mat files for the analysis of the residence time',...
        'MultiSelect','on');
end
if ~iscell(FileNames)
    if FileNames == 0
        FileNames = [];
    end
end
% Set the parameters for the analysis of the residence time;

if ~isempty(FileNames)
    if isfield(handles.Parameters, 'Used') && isfield(handles.Parameters.Used, 'Bound');
        
        defPar(1) = handles.Parameters.Used.Bound.ThreshL;
        defPar(2) = handles.Parameters.Used.Bound.ThreshH;
        defPar(3) = handles.Parameters.Used.Bound.minBoundFrames;
        defPar(4) = handles.Parameters.Used.Acquisition.frameTime;
        defPar(5) = 0;
        
    else
        defPar(1) = handles.Parameters.Bound.ThreshL;
        defPar(2) = handles.Parameters.Bound.ThreshH;
        defPar(3) = handles.Parameters.Bound.minBoundFrames;
        defPar(4) = handles.Parameters.Acquisition.frameTime;
        defPar(5) = 0;
        
    end;
    
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
        num2str(defPar(3)), '0', num2str(defPar(4)), num2str(defPar(5))};
    prompt = {'Maximum Jump between consecutive frames',...
        'Maximum end to end distance',...
        'Minimum length of bound tracks',...
        'Points in the histogram after log sampling (0 = No Log Sampling)',...
        'Frame Time',...
        'Maximum Frames to Analyze (0 = all)'};
    dlgtitle = 'Parameters for the Survival Time distribution';
    AnalysisParam = inputdlg(prompt,dlgtitle, 1, defaults);
    
    %     end
    
    % Send parameters and file names to the residence time GUI
    AnalyzeResTimes_GUI(AnalysisFlag,PathName,FileNames, AnalysisParam);
end







% --------------------------------------------------------------------
function JDH_Merge_Callback(hObject, eventdata, handles)
% hObject    handle to JDH_Merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Select multiple Files to open;

if isfield(handles, 'Data')
    [FileNames,PathName,FilterIndex] = uigetfile('.mat',...
        'Select Mat files for the analysis of the residence time',...
        [handles.Data.pathName, handles.Data.fileName],'MultiSelect','on');
else
    [FileNames,PathName,FilterIndex] = uigetfile('.mat',...
        'Select Mat files for the analysis of the residence time',...
        'MultiSelect','on');
end
if ~iscell(FileNames)
    if FileNames == 0
        FileNames = [];
    end
end
if ~isempty(FileNames)
    if isfield(handles.Parameters, 'Used') && isfield(handles.Parameters.Used, 'Bound');
        
        defPar(1) = handles.Parameters.Used.Analysis.JHistBinSize;
        defPar(2) = handles.Parameters.Used.Analysis.JHistMaxFrameN;
        defPar(3) = handles.Parameters.Used.Tracking(5).*...
            handles.Parameters.Used.Acquisition.pixelSize;
        defPar(4) = handles.Parameters.Used.Acquisition.frameTime;
        
        
    else
        defPar(1) = handles.Parameters.Analysis.JHistBinSize;
        defPar(2) = handles.Parameters.Analysis.JHistMaxFrameN;
        defPar(3) = 12*handles.Parameters.Acquisition.pixelSize;
        defPar(4) = handles.Parameters.Acquisition.frameTime;
        
    end;
    
    % Get a dialog input to set the parameters;
    
    defaults = {num2str(defPar(1)), num2str(defPar(2)), ...
        num2str(defPar(3)), num2str(defPar(4))};
    prompt = {'Size of the histogram bin [um]',...
        'Number of timepoints for the histogram',...
        'Maximum distance [um]',...
        'Frame Time'};
    
    dlgtitle = 'Parameters for histogram of displacements';
    AnalysisParam = inputdlg(prompt,dlgtitle, 1, defaults);
    AnalyzeHistD_GUI([],PathName,FileNames, AnalysisParam);
end


% --- Executes on selection change in ROIList.
function ROIList_Callback(hObject, eventdata, handles)
% hObject    handle to ROIList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ROIList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ROIList
curSel = get(hObject,'Value');

plotROI(handles.Process.ROIpos,handles.axes1,curSel);
if ~isempty(handles.Process.ROIClass) && ~isempty(handles.Process.ROIClass{1,:}) && ~isempty(handles.Process.ROIClass{curSel,:})
    if iscell(handles.Process.ROIClass{1,:})
        set(handles.ROIClassCur_text,'String',handles.Process.ROIClass{curSel,:});
    else
        set(handles.ROIClassCur_text,'String','Not Defined');
    end
else
    set(handles.ROIClassCur_text,'String','Not Defined');
end


% --- Executes during object creation, after setting all properties.
function ROIList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROIList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in roiRemovePush.
function roiRemovePush_Callback(hObject, eventdata, handles)
% hObject    handle to roiRemovePush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

answer = questdlg('Do you want to delete the selected ROI?','Delete ROI','Yes','No','No');

switch answer
    case 'Yes'
        ROI_selected = get(handles.ROIList,'Value');
        %Check if tracking data exists
        if isfield(handles,'Tracking')
            ROI_nest_mat = roiNest(handles.Process.ROIimage);
            if max(ROI_nest_mat(:,ROI_selected)) == 1 && (isfield(handles.Tracking,'Centroids') || isfield(handles.Tracking,'Particles'))
                reassignPartAns = questdlg('ROI contains Particles, Should they be reassigned to the enclosing ROI or Deleted?', ...
                    'Particles in ROI','Reassign','Delete','Cancel','Reassign');
                switch reassignPartAns
                    case 'Reassign'
                        newROI_ind = find(ROI_nest_mat(:,ROI_selected) == 1,1);
                        if isfield(handles.Tracking,'Centroids')
                            Centroids = handles.Tracking.Centroids;
                            Centroids(Centroids(:,7) == ROI_selected, 7) = newROI_ind;
                            Centroids(Centroids(:,7) > ROI_selected,7) = Centroids(Centroids(:,7) > ROI_selected,7) - 1;
                            Centroids = sortrows(Centroids,[6 7]);
                            handles.Tracking.Centroids = Centroids;
                        end
                        
                        if isfield(handles.Tracking,'Particles')
                            Particles = handles.Tracking.Particles;
                            Particles(Particles(:,13) == ROI_selected,13) = newROI_ind;
                            Particles(Particles(:,13) > ROI_selected,13) = Particles(Particles(:,13) > ROI_selected,13) - 1;
                            Particles = sortrows(Particles,[6 13]);
                            handles.Tracking.Particles = Particles;
                        end
                        
                        if isfield(handles.Tracking,'Tracks')
                            Tracks = handles.Tracking.Tracks;
                            Tracks(Tracks(:,5) == ROI_selected,5) = newROI_ind;
                            Tracks(Tracks(:,5) > ROI_selected, 5) = Tracks(Tracks(:,5) > ROI_selected, 5) - 1;
                            handles.Tracking.Tracks = Tracks;
                            
                        end
                        if isfield(handles.Tracking,'CheckTracks')
                            Tracks = handles.Tracking.CheckTracks;
                            Tracks(Tracks(:,5) == ROI_selected,5) = newROI_ind;
                            Tracks(Tracks(:,5) > ROI_selected, 5) = Tracks(Tracks(:,5) > ROI_selected, 5) - 1;
                            handles.Tracking.CheckTracks = Tracks;
                            
                        end
                        if isfield(handles,'PreAnalysis')
                            if isfield(handles.PreAnalysis,'Tracks_um')
                                Tracks = handles.PreAnalysis.Tracks_um;
                                Tracks(Tracks(:,5) == ROI_selected,5) = newROI_ind;
                                Tracks(Tracks(:,5) > ROI_selected,5) = Tracks(Tracks(:,5) > ROI_selected,5) - 1;
                                handles.PreAnalysis.Tracks_um = Tracks;
                            end
                            if isfield(handles.PreAnalysis,'NParticles')
                                if size(handles.PreAnalysis.NParticles,2) >= ROI_selected+1
                                    NParticles = handles.PreAnalysis.NParticles;
                                    NParticles(:,newROI_ind+1) = NParticles(:,newROI_ind+1) + NParticles(:,ROI_selected+1);
                                    NParticles(:,ROI_selected+1) = [];
                                    handles.PreAnalysis.NParticles = NParticles;
                                end
                            end
                        end
                        
                        if isfield(handles.Tracking,'CheckParticles');
                            Particles = handles.Tracking.CheckParticles;
                            if size(Particles,2) > 8
                                ROI_col = 13;
                            else
                                ROI_col = 7;
                            end
                            Particles(Particles(:,ROI_col) == ROI_selected,13) = newROI_ind;
                            Particles(Particles(:,ROI_col) > ROI_selected,13) = Particles(Particles(:,ROI_col) > ROI_selected,13) - 1;
                            Particles = sortrows(Particles,[6,ROI_col]);
                            handles.Tracking.CheckParticles = Particles;
                        end
                    case 'Delete'
                        if isfield(handles.Tracking,'Centroids')
                            Centroids = handles.Tracking.Centroids;
                            Centroids(Centroids(:,7) == ROI_selected, :) = [];
                            Centroids(Centroids(:,7) > ROI_selected,7) = Centroids(Centroids(:,7) > ROI_selected,7) - 1;
                            Centroids = sortrows(Centroids,[6 7]);
                            handles.Tracking.Centroids = Centroids;
                        end
                        
                        if isfield(handles.Tracking,'Particles')
                            Particles = handles.Tracking.Particles;
                            Particles(Particles(:,13) == ROI_selected,:) = [];
                            Particles(Particles(:,13) > ROI_selected,13) = Particles(Particles(:,13) > ROI_selected,13) - 1;
                            Particles = sortrows(Particles,[6 13]);
                            handles.Tracking.Particles = Particles;
                        end
                        
                        if isfield(handles.Tracking,'Tracks')
                            Tracks = handles.Tracking.Tracks;
                            Tracks(Tracks(:,5) == ROI_selected,:) = [];
                            Tracks(Tracks(:,5) > ROI_selected, 5) = Tracks(Tracks(:,5) > ROI_selected, 5) - 1;
                            handles.Tracking.Tracks = Tracks;
                            
                        end
                        if isfield(handles.Tracking,'CheckTracks')
                            Tracks = handles.Tracking.CheckTracks;
                            Tracks(Tracks(:,5) == ROI_selected,:) = [];
                            Tracks(Tracks(:,5) > ROI_selected, 5) = Tracks(Tracks(:,5) > ROI_selected, 5) - 1;
                            handles.Tracking.CheckTracks = Tracks;
                            
                        end
                        if isfield(handles,'PreAnalysis')
                            if isfield(handles.PreAnalysis,'Tracks_um')
                                Tracks = handles.PreAnalysis.Tracks_um;
                                Tracks(Tracks(:,5) == ROI_selected,:) = [];
                                Tracks(Tracks(:,5) > ROI_selected,5) = Tracks(Tracks(:,5) > ROI_selected,5) - 1;
                                handles.PreAnalysis.Tracks_um = Tracks;
                            end
                            if isfield(handles.PreAnalysis,'NParticles')
                                if size(handles.PreAnalysis.NParticles,2) >= ROI_selected+1
                                    NParticles = handles.PreAnalysis.NParticles;
                                    NParticles(:,ROI_selected+1) = [];
                                    handles.PreAnalysis.NParticles = NParticles;
                                end
                                
                            end
                        end
                        
                        if isfield(handles.Tracking,'CheckParticles');
                            Particles = handles.Tracking.CheckParticles;
                            if size(Particles,2) > 8
                                ROI_col = 13;
                            else
                                ROI_col = 7;
                            end
                            Particles(Particles(:,ROI_col) == ROI_selected,:) = [];
                            Particles(Particles(:,ROI_col) > ROI_selected,13) = Particles(Particles(:,ROI_col) > ROI_selected,13) - 1;
                            Particles = sortrows(Particles,[6,ROI_col]);
                            handles.Tracking.CheckParticles = Particles;
                        end
                    case 'Cancel'
                        return
                        
                end
            elseif max(ROI_nest_mat(:,ROI_selected)) == 0 && (isfield(handles.Tracking,'Centroids') || isfield(handles.Tracking,'Particles'))
                deletePartAns = questdlg('Deleting the ROI will remove Particle data','Delete Particles','OK','Cancel','Cancel');
                switch deletePartAns
                    case 'OK'
                        if length(handles.Process.ROIpos) == 1
                            if isfield(handles.Tracking,'Centroids')
                                handles.Tracking = rmfield(handles.Tracking,'Centroids');
                            end
                            if isfield(handles.Tracking,'Particles')
                                handles.Tracking = rmfield(handles.Tracking,'Particles');
                            end
                            if isfield(handles.Tracking,'CheckParticles')
                                handles.Tracking = rmfield(handles.Tracking,'CheckParticles');
                            end
                            if isfield(handles.Tracking,'Tracks')
                                handles.Tracking = rmfield(handles.Tracking,'Tracks');
                            end
                            if isfield(handles,'PreAnalysis')
                                handles = rmfield(handles,'PreAnalysis');
                            end
                        else
                            if isfield(handles.Tracking,'Centroids')
                                Centroids = handles.Tracking.Centroids;
                                Centroids(Centroids(:,7) == ROI_selected, :) = [];
                                Centroids(Centroids(:,7) > ROI_selected,7) = Centroids(Centroids(:,7) > ROI_selected,7) - 1;
                                Centroids = sortrows(Centroids,[6 7]);
                                handles.Tracking.Centroids = Centroids;
                            end
                            
                            if isfield(handles.Tracking,'Particles')
                                Particles = handles.Tracking.Particles;
                                Particles(Particles(:,13) == ROI_selected,:) = [];
                                Particles(Particles(:,13) > ROI_selected,13) = Particles(Particles(:,13) > ROI_selected,13) - 1;
                                Particles = sortrows(Particles,[6 13]);
                                handles.Tracking.Particles = Particles;
                            end
                            if isfield(handles.Tracking,'Tracks')
                                Tracks = handles.Tracking.Tracks;
                                Tracks(Tracks(:,5) == ROI_selected,:) = [];
                                Tracks(Tracks(:,5) > ROI_selected, 5) = Tracks(Tracks(:,5) > ROI_selected, 5) - 1;
                                handles.Tracking.Tracks = Tracks;
                                
                            end
                            if isfield(handles.Tracking,'CheckTracks')
                                Tracks = handles.Tracking.CheckTracks;
                                Tracks(Tracks(:,5) == ROI_selected,:) = [];
                                Tracks(Tracks(:,5) > ROI_selected, 5) = Tracks(Tracks(:,5) > ROI_selected, 5) - 1;
                                handles.Tracking.CheckTracks = Tracks;
                                
                            end
                            if isfield(handles,'PreAnalysis')
                                if isfield(handles.PreAnalysis,'Tracks_um')
                                    Tracks = handles.PreAnalysis.Tracks_um;
                                    Tracks(Tracks(:,5) == ROI_selected,:) = [];
                                    Tracks(Tracks(:,5) > ROI_selected,5) = Tracks(Tracks(:,5) > ROI_selected,5) - 1;
                                    handles.PreAnalysis.Tracks_um = Tracks;
                                end
                                if isfield(handles.PreAnalysis,'NParticles')
                                    NParticles = handles.PreAnalysis.NParticles;
                                    NParticles(:,ROI_selected+1) = [];
                                    handles.PreAnalysis.NParticles = NParticles;
                                    
                                end
                            end
                            
                            if isfield(handles.Tracking,'CheckParticles');
                                Particles = handles.Tracking.CheckParticles;
                                if size(Particles,2) > 8
                                    ROI_col = 13;
                                else
                                    ROI_col = 7;
                                end
                                Particles(Particles(:,ROI_col) == ROI_selected,:) = [];
                                Particles(Particles(:,ROI_col) > ROI_selected,13) = Particles(Particles(:,ROI_col) > ROI_selected,13) - 1;
                                Particles = sortrows(Particles,[6,ROI_col]);
                                handles.Tracking.CheckParticles = Particles;
                            end
                            
                        end
                    case 'Cancel'
                        return
                        
                end
            end
        end
        %Delete the ROI from the listbox
        ROI_list = get(handles.ROIList,'String');
        if length(ROI_list) > 1
            
            
            ROI_listb = cell(length(ROI_list)-1,1);
            ind = 1;
            for i = 1:length(ROI_list)
                if i ~= ROI_selected
                    ROI_listb{ind,:} = ROI_list{i,:};
                    ind = ind +1;
                end
            end
            set(handles.ROIList,'String',ROI_listb);
            set(handles.ROIList,'Value',1);
            
            %Remove the ROI data from the GUI data
            if isfield(handles.Process,'ROIpos')
                ind = 1;
                ROIpos_tmp = cell(length(handles.Process.ROIpos)-1,1);
                ROIimage_tmp = cell(length(handles.Process.ROIpos)-1,1);
                ROIlabel_tmp = cell(length(handles.Process.ROIpos)-1,1);
                for i = 1:length(handles.Process.ROIpos)
                    if i ~= ROI_selected
                        ROIpos_tmp{ind,:} = handles.Process.ROIpos{i,:};
                        ROIimage_tmp{ind,:} = handles.Process.ROIimage{i,:};
                        ROIlabel_tmp{ind,:} = handles.Process.ROIlabel{i,:};
                        ind = ind + 1;
                    end
                end
                handles.Process.ROIpos = ROIpos_tmp;
                handles.Process.ROIimage = ROIimage_tmp;
                handles.Process.ROIlabel = ROIlabel_tmp;
            end
            %Remove the class
            if ~isempty(handles.Process.ROIClass)
                if ~isempty(handles.Process.ROIClass{1,:}) && iscell(handles.Process.ROIClass{1,:})
                    ind = 1;
                    ROIClass_tmp = cell(length(handles.Process.ROIClass)-1,1);
                    for i = 1:length(handles.Process.ROIClass)
                        if i ~= ROI_selected
                            ROIClass_tmp{ind,:} = handles.Process.ROIClass{i,:};
                            
                            ind = ind+1;
                        else
                            ROIClass_rem = handles.Process.ROIClass{i,:};
                        end
                    end
                    handles.Process.ROIClass = ROIClass_tmp;
                    remClass = 1;
                    for i = 1:length(handles.Process.ROIClass)
                        if strcmp(ROIClass_rem,handles.Process.ROIClass{i,:});
                            remClass = 0;
                            break
                        end
                    end
                    if remClass == 1
                        
                        for i = 1:length(handles.Process.AllROIClasses)
                            if strcmp(ROIClass_rem,handles.Process.AllROIClasses{i,:})
                                handles.Process.AllROIClasses{i,:} = [];
                                break
                            end
                        end
                    end
                    
                    
                end
            end
            
            %Replot the image, remaining ROIs, and particles
            set(handles.axes1,'NextPlot','replacechildren');
            imagesc(handles.Current.Image, handles.Current.clims);
            %             axis image;
            plotROI(handles.Process.ROIpos);
            if get(handles.showParticles,'Value')
                if handles.isFitPSF
                    plotParticle(handles.Tracking.Particles, ...
                        handles.Current.Ix, handles.isFitPSF);
                else
                    plotParticle(handles.Tracking.Centroids, ...
                        handles.Current.Ix, handles.isFitPSF);
                end
            end
            
        else
            %if this is the last remaining ROI, make it appear as it did at
            %the beginning
            set(handles.ROIList,'String','');
            set(handles.ROIList,'Value',1);
            handles.Process = rmfield(handles.Process,{'ROIpos','ROIimage','ROIlabel'});
            handles.Process.ROIClass = [];
            if isfield(handles.Process,'AllROIClasses');
                handles.Process = rmfield(handles.Process,'AllROIClasses');
            end
            set(handles.roiRemovePush,'Enable','off');
            set(handles.roiNamePush,'Enable','off');
            set(handles.CopyROI,'Enable','off');
            %Replot the image & remaining ROIs
            set(handles.axes1,'NextPlot','replacechildren');
            imagesc(handles.Current.Image, handles.Current.clims);
            %             axis image;
        end
        
        
        % Update handles structure
        guidata(hObject, handles);
        
    case 'No'
        return
end

guidata(hObject,handles);



% --- Executes on button press in roiNamePush.
function roiNamePush_Callback(hObject, eventdata, handles)
% hObject    handle to roiNamePush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ROI_list = get(handles.ROIList,'String');
ROIselected = get(handles.ROIList,'Value');

def_name = {ROI_list{ROIselected,:}};

ROI_label = inputdlg('Specify a new label for the selected ROI:','ROI label',1,def_name);

ROI_list{ROIselected,:} = ROI_label{1,:};

set(handles.ROIList,'String',ROI_list);

handles.Process.ROIlabel = ROI_list;

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function TwoColorSMC_Callback(hObject, eventdata, handles)
% hObject    handle to TwoColorSMC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'Tracking')
    if isfield(handles.Tracking,'Particles')
        if ~isempty(handles.Tracking.Particles)
            UseCurrentAns = questdlg('Use the current analysis as 1 of the channels?', 'Export Data','Yes','No','Cancel','Yes');
            switch UseCurrentAns
                case 'Yes'
                    TwoColorSMcoloc_GUI(handles.Data,handles.Process,handles.Tracking,handles.PreAnalysis,handles.Parameters);
                    
                case 'No'
                    TwoColorSMcoloc_GUI();
                case 'Cancel'
                    return
            end
        else
            TwoColorSMcoloc_GUI();
        end
    else
        TwoColorSMcoloc_GUI();
    end
else
    TwoColorSMcoloc_GUI();
end




% --------------------------------------------------------------------
function ImportSOS_Callback(hObject, eventdata, handles)
% hObject    handle to ImportSOS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Check that an image is loaded and ROI is selected
if ~isfield(handles,'Data')
    msgbox('Load in Image data first','No image data');
    return
else
    if ~isfield(handles.Data,'imageStack')
        msgbox('Load in Image data first','No image data');
        return
    else
        if isempty(handles.Data.imageStack)
            msgbox('Load in Image data first','No image data');
            return
        end
    end
    
end

if ~isfield(handles,'Process')
    msgbox('Specify ROI(s) first','No ROI data');
    return
else
    if  ~isfield(handles.Process,'ROIpos')
        msgbox('Specify ROI(s) first','No ROI data');
        return
    else
        if isempty(handles.Process.ROIpos)
            msgbox('Specify ROI(s) first','No ROI data');
            return
        end
    end
end

if ~isempty(handles.Tracking);
    answer = questdlg('Tracking data already available. Do you want to overwrite it?',...
        'Overwrite?','Keep', 'Overwrite', 'Keep');
else
    answer = 'Overwrite';
end

switch answer
    case 'Overwrite'
        
        set(handles.showTracks, 'Value', 0);
        set(handles.showTracks, 'Enable', 'off');
        delete(findobj(gca,'Color','g'));
        set(handles.TrackPopUp, 'Value',1);
        set(handles.TrackPopUp, 'Enable', 'off');
        
        % Empty vectors with Tracking and Analysis data;
        handles.Tracking = [];
        handles.PreAnalysis = [];
        handles.Analysis = [];
        
        
        
        
        
        if isfield(handles.Data,'pathName')
            defDir = handles.Data.pathName;
        else
            defDir = pwd;
        end
        
        SOSfolder = uigetdir(defDir,'Select folder containing SOS detection & tracking results');
        
        if SOSfolder ~= 0
            folderContents = dir(SOSfolder);
            foundDetections = 0;
            foundTracking = 0;
            foundParamters = 0;
            %Go through the folder contents searching for the appropriate files
            for i = 1:length(folderContents)
                if strcmp(folderContents(i).name,'detections.txt')
                    foundDetections = i;
                end
                if strcmp(folderContents(i).name,'tracks.simple.filtered.txt')
                    foundTracking = i;
                end
                if strcmp(folderContents(i).name,'parameters.txt');
                    foundParameters = i;
                end
                %If they have both been found, we can stop the search
                if foundTracking > 0 && foundDetections > 0 && foundParameters > 0
                    break
                end
            end
            %if we don't find any of these files, we can't continue
            if foundDetections == 0
                msgbox('Could not locate detections.txt in selected folder', 'No Particle Data');
                return
            end
            
            if foundTracking == 0
                msgbox('Could not locate tracking.simple.filtered.txt in selected folder', 'No Tracking Data');
                return
            end
            
            if foundParameters == 0
                msgbox('Could not locate parameters.txt in selected folder', 'No Tracking Data');
                return
            end
            %specify data size & format
            sizeP = [10 Inf];
            sizeT = [8 Inf];
            
            formatSpec = '%f';
            
            % read in Particle Data
            fidP = fopen([SOSfolder,filesep,'detections.txt'],'r');
            SOSParticles = fscanf(fidP,formatSpec,sizeP);
            SOSParticles = SOSParticles';
            fclose(fidP);
            
            %read in Tracking Data
            fidT = fopen([SOSfolder,filesep,'tracks.simple.filtered.txt'],'r');
            SOSTracks = fscanf(fidT,formatSpec,sizeT);
            SOSTracks = SOSTracks';
            fclose(fidT);
            
            %read in Parameters
            fid = fopen([SOSfolder,filesep,'parameters.txt'],'r');
            tline1 = fgetl(fid); %first line of file is text, which we don't need
            SOSParams = fscanf(fid,formatSpec);
            fclose(fid);
            
            handles.Parameters.Used.Tracking(1) = SOSParams(7); %loBP (~s_min)
            handles.Parameters.Used.Tracking(2) = SOSParams(8); %hiBP (~s_max)
            handles.Parameters.Used.Tracking(3) = SOSParams(5); %Threshold (I_min)
            handles.Parameters.Used.Tracking(4) = SOSParams(10); %windowSize(masksize)
            handles.Parameters.Used.Tracking(5) = 3*SOSParams(3); %maxJump (sigma_diffusion)
            handles.Parameters.Used.Tracking(6) = 1; %fillGaps (always 1)
            handles.Parameters.Used.Tracking(7) = SOSParams(15); %shortestTrack (deltaT)
            set(handles.thresholdEdit,'String',num2str(SOSParams(5)));
            %initialize a vector to put the particle data from SOS into the same
            %format expected here
            Particles = zeros(size(SOSParticles,1),13);
            
            Particles(:,1) = SOSParticles(:,1)+1; %x position
            Particles(:,2) = SOSParticles(:,2)+1; %y position
            Particles(:,3) = SOSParticles(:,4); %Integrated Intensity
            %     Particles(:,4) = SOSParticles(:,2); %radius of gyration (don't think its included)
            %     Particles(:,5) = SOSParticles(:,6); %center intensity (guess)
            Particles(:,6) = SOSParticles(:,3) + 1;
            %     Particles(:,7) = SOSParticles(:,9); %amplitude of Gaussian fit (guess)
            Particles(:,8) = SOSParticles(:,5); %sigma
            %     Particles(:,9) = SOSParticles(:,4); %Background (don't think its included)
            Particles(:,10) = SOSParticles(:,1)+1; %X position (Gauss fit)
            Particles(:,11) = SOSParticles(:,2)+1; %Y position (Gauss fit)
            Particles(:,12) = SOSParticles(:,6); %SSR of fit (don't think it's included)
            
            
            
            
            %initialize a vector to put the tracking data from SOS into the same
            %format expected here
            Tracks = zeros(size(SOSTracks,1),5);
            
            Tracks(:,1) = SOSTracks(:,2)+1; % x position
            Tracks(:,2) = SOSTracks(:,3)+1; % y position
            Tracks(:,3) = SOSTracks(:,1) + 1; % frame number
            Tracks(:,4) = SOSTracks(:,4) + 1; % Particle Number
            
            %ensure that there are no gaps in the tracks
            Tracks_all = [];
            for i = 1:max(Tracks(:,4))
                if i == 6
                    boa = 10;
                end
                
                curTrack = Tracks(Tracks(:,4) == i,:);
                curTrack2 = [];
                if ~isempty(curTrack)
                    minFrame = min(curTrack(:,3));
                    maxFrame = max(curTrack(:,3));
                    for j = minFrame:maxFrame
                        if isempty(curTrack(curTrack(:,3) == j,:))
                            lastPt = j-1;
                            lastJmp = 1;
                            while isempty(curTrack(curTrack(:,3) == lastPt,:))
                                lastPt = lastPt - 1;
                                lastJmp = lastJmp + 1;
                            end
                            nextPt = j+1;
                            nextJmp = 1;
                            while isempty(curTrack(curTrack(:,3) == nextPt,:))
                                nextPt = nextPt + 1;
                                nextJmp = nextJmp + 1;
                            end
                            AddPt(1,3) = j;
                            AddPt(1,4) = i;
                            AddPt(1,5) = 0;
                            lastX = curTrack(curTrack(:,3) == lastPt,1);
                            lastY = curTrack(curTrack(:,3) == lastPt,2);
                            nextX = curTrack(curTrack(:,3) == nextPt,1);
                            nextY = curTrack(curTrack(:,3) == nextPt,2);
                            AddPt(1,1) = ((lastX*nextJmp) + (nextX*lastJmp))/(nextJmp+lastJmp);
                            AddPt(1,2) = ((lastY*nextJmp) + (nextY*lastJmp))/(nextJmp+lastJmp);
                            curTrack2 = [curTrack2; AddPt];
                            Particles2Add = [AddPt(1,1:2), 0, 0, 0, AddPt(1,3), 0,0,0,AddPt(1,1:2),0,0];
                            Particles = [Particles;Particles2Add];
                        else
                            curTrack2 = [curTrack2; curTrack(curTrack(:,3) == j,:)];
                        end
                    end
                end
                Tracks_all = [Tracks_all; curTrack2];
            end
            Particles2 = InsideROIcheck2(Particles,handles.Process.ROIimage);
            
            handles.Tracking.Particles = Particles2;
            Tracks2 = InsideROIcheck2(Tracks_all,handles.Process.ROIimage);
            
            handles.Tracking.Tracks = Tracks2;
            handles.Current.Tracks = handles.Tracking.Tracks;
            set(handles.showParticles,'Enable','on');
            set(handles.showTracks,'Enable','on');
            handles.isFitPSF = 1;
            
            set(handles.checkButton,'Enable','on');
            set(handles.TrackDimButton,'Enable','on');
            set(handles.Track_Preprocess,'Enable','on');
            set(handles.HistD_Analysis,'Enable','on');
            set(handles.Bound_Analysis,'Enable','on');
            set(handles.Export_Params,'Enable','on');
        end
end

guidata(hObject,handles);


% --- Executes on button press in LoadROIfromFile.
function LoadROIfromFile_Callback(hObject, eventdata, handles)
% hObject    handle to LoadROIfromFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName] = ...
    uigetfile('.mat','Select Matlab File with ROI data');

if FileName ~= 0
    set(handles.StatusText,'String','Reading MAT file to extract ROIs...');
    drawnow;
    ListVar = whos('-file',[PathName,FileName]);
    if ismember('Results',{ListVar.name});
        IN = load([PathName, FileName]);
        if isfield(IN.Results.Process,'ROIpos')
            if ~isempty(IN.Results.Process.ROIpos)
                handles.Process.ROIpos = IN.Results.Process.ROIpos;
                handles.Process.ROIimage = IN.Results.Process.ROIimage;
                handles.Process.ROIlabel = IN.Results.Process.ROIlabel;
                %Add classes if they exist
                if isfield(IN.Results.Process,'ROIClass');
                    handles.Process.ROIClass = IN.Results.Process.ROIClass;
                    if ~isempty(IN.Results.Process.ROIClass)
                        if ~isempty(IN.Results.Process.ROIClass{1,:}) && iscell(IN.Results.Process.ROIClass{1,:})
                            handles.Process.AllROIClasses = IN.Results.Process.AllROIClasses;
                        end
                    end
                end
                set(handles.roiNamePush,'Enable','on');
                set(handles.roiRemovePush,'Enable','on');
                nROIs = length(handles.Process.ROIpos);
                for i = 1:nROIs
                    ROI_list{i,:} = handles.Process.ROIlabel{i,:};
                    plotROI(handles.Process.ROIpos);
                end
                set(handles.ROIList,'String',ROI_list);
                
            end
        end
    end
    
    set(handles.StatusText,'String','Reading MAT file to extract ROIs...Done');
    drawnow;
    
end

guidata(hObject,handles);


% --------------------------------------------------------------------
function uTrackImport_Callback(hObject, eventdata, handles)
% hObject    handle to uTrackImport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Check that an image is loaded and ROI is selected
if ~isfield(handles,'Data')
    msgbox('Load in Image data first','No image data');
    return
else
    if ~isfield(handles.Data,'imageStack')
        msgbox('Load in Image data first','No image data');
        return
    else
        if isempty(handles.Data.imageStack)
            msgbox('Load in Image data first','No image data');
            return
        end
    end
    
end

if ~isfield(handles,'Process')
    msgbox('Specify ROI(s) first','No ROI data');
    return
else
    if  ~isfield(handles.Process,'ROIpos')
        msgbox('Specify ROI(s) first','No ROI data');
        return
    else
        if isempty(handles.Process.ROIpos)
            msgbox('Specify ROI(s) first','No ROI data');
            return
        end
    end
end

if ~isempty(handles.Tracking);
    answer = questdlg('Tracking data already available. Do you want to overwrite it?',...
        'Overwrite?','Keep', 'Overwrite', 'Keep');
else
    answer = 'Overwrite';
end

switch answer
    case 'Overwrite'
        set(handles.showTracks, 'Value', 0);
        set(handles.showTracks, 'Enable', 'off');
        delete(findobj(gca,'Color','g'));
        set(handles.TrackPopUp, 'Value',1);
        set(handles.TrackPopUp, 'Enable', 'off');
        
        % Empty vectors with Tracking and Analysis data;
        handles.Tracking = [];
        handles.PreAnalysis = [];
        handles.Analysis = [];
        
        
        
        
        
        if isfield(handles.Data,'pathName')
            defDir = handles.Data.pathName;
        else
            defDir = pwd;
        end
        
        uTfolder = uigetdir(defDir,'Select folder containing uTrack detection & tracking results');
        
        if uTfolder ~= 0
            
            ParticleFile = [uTfolder, filesep, 'GaussianMixtureModels',filesep,'Channel_1_detection_result.mat'];
            TrackFile = [uTfolder,filesep,'tracks',filesep,'Channel_1_tracking_result.mat'];
            
            PartIN = load(ParticleFile);
            TrackIN = load(TrackFile);
            
            PartStruct = PartIN.movieInfo;
            TrackStruct = TrackIN.tracksFinal;
            handles.Parameters.Used.Tracking(1) = 0; %loBP (~s_min)
            handles.Parameters.Used.Tracking(2) = 100; %hiBP (~s_max)
            handles.Parameters.Used.Tracking(3) = 0; %Threshold (I_min)
            handles.Parameters.Used.Tracking(4) = 0; %windowSize(masksize)
            handles.Parameters.Used.Tracking(5) = 4; %maxJump (sigma_diffusion)
            Particles = [];
            for i = 1:size(PartStruct,1)
                Particles_tmp = zeros(size(PartStruct(i).xCoord,1),13);
                for j = 1:size(PartStruct(i).xCoord,1)
                    Particles_tmp(j,1) = PartStruct(i).xCoord(j,1);
                    Particles_tmp(j,2) = PartStruct(i).yCoord(j,1);
                    Particles_tmp(j,3) = PartStruct(i).amp(j,1);
                    Particles_tmp(j,6) = i;
                    Particles_tmp(j,8) = PartStruct(i).xCoord(j,2);
                    Particles_tmp(j,10) = PartStruct(i).xCoord(j,1);
                    Particles_tmp(j,11) = PartStruct(i).yCoord(j,1);
                end
                Particles = [Particles; Particles_tmp];
                
            end
            
            
            Tracks = [];
            for i = 1:size(TrackStruct,1)
                nPts = size(TrackStruct(i).tracksFeatIndxCG,2);
                TrackPtData = reshape(TrackStruct(i).tracksCoordAmpCG,8,nPts);
                TrackPtData = TrackPtData';
                tVec = (TrackStruct(i).seqOfEvents(TrackStruct(i).seqOfEvents(:,2) == 1,1):TrackStruct(i).seqOfEvents(TrackStruct(i).seqOfEvents(:,2) == 2,1))';
                Tracks_tmp = zeros(nPts,5);
                
                Tracks_tmp(:,1) = TrackPtData(:,1); % x position
                Tracks_tmp(:,2) = TrackPtData(:,2); % y position
                Tracks_tmp(:,3) = tVec; % frame number
                Tracks_tmp(:,4) = i; % Particle Number
                Tracks = [Tracks; Tracks_tmp];
            end
            Particles2 = InsideROIcheck2(Particles,handles.Process.ROIimage);
            
            handles.Tracking.Particles = Particles2;
            Tracks2 = InsideROIcheck2(Tracks,handles.Process.ROIimage);
            
            handles.Tracking.Tracks = Tracks2;
            handles.Current.Tracks = handles.Tracking.Tracks;
            set(handles.showParticles,'Enable','on');
            set(handles.showTracks,'Enable','on');
            handles.isFitPSF = 1;
            
            set(handles.checkButton,'Enable','on');
            set(handles.TrackDimButton,'Enable','on');
            set(handles.Track_Preprocess,'Enable','on');
            set(handles.HistD_Analysis,'Enable','on');
            set(handles.Bound_Analysis,'Enable','on');
            set(handles.Export_Params,'Enable','on');
            
            guidata(hObject,handles);
        end
end





% --------------------------------------------------------------------
function ResTime_Merge2_Callback(hObject, eventdata, handles)
% hObject    handle to ResTime_Merge2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Select multiple Files to open;

if isfield(handles, 'Data')
    [FileNames1,PathName1,FilterIndex] = uigetfile('.mat',...
        'Select Mat files for the analysis of the residence time',...
        [handles.Data.pathName, handles.Data.fileName],'MultiSelect','on');
else
    [FileNames1,PathName1,FilterIndex] = uigetfile('.mat',...
        'Select Mat files for the analysis of the residence time',...
        'MultiSelect','on');
end
if ~iscell(FileNames1)
    if FileNames1 == 0
        FileNames1 = [];
    end
end
% Set the parameters for the analysis of the residence time;
if ~isempty(FileNames1)
    if isfield(handles.Parameters, 'Used') && isfield(handles.Parameters.Used, 'Bound');
        
        defPar(1) = handles.Parameters.Used.Bound.ThreshL;
        defPar(2) = handles.Parameters.Used.Bound.ThreshH;
        defPar(3) = handles.Parameters.Used.Bound.minBoundFrames;
        defPar(4) = handles.Parameters.Used.Acquisition.frameTime;
        
    else
        defPar(1) = handles.Parameters.Bound.ThreshL;
        defPar(2) = handles.Parameters.Bound.ThreshH;
        defPar(3) = handles.Parameters.Bound.minBoundFrames;
        defPar(4) = handles.Parameters.Acquisition.frameTime;
        
    end;
    
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
    if isfield(handles, 'Data')
        [FileNames2,PathName2,FilterIndex] = uigetfile('.mat',...
            'Select Mat files for the analysis of the residence time',...
            [handles.Data.pathName, handles.Data.fileName],'MultiSelect','on');
    else
        [FileNames2,PathName2,FilterIndex] = uigetfile('.mat',...
            'Select Mat files for the analysis of the residence time',...
            'MultiSelect','on');
    end
    if ~iscell(FileNames2)
        if FileNames2 == 0
            FileNames2 = [];
        end
    end
    if ~isempty(FileNames2)
        defaults = {num2str(defPar(4))};
        prompt = {'Frame Time'};
        dlgtitle = 'Parameters for the Second set of Data';
        AnalysisParam2 = inputdlg(prompt,dlgtitle, 1, defaults);
        AnalysisParam{6} = AnalysisParam2{1};
        % Send parameters and file names to the residence time GUI
        AnalyzeResTimesDiffRates_GUI(AnalysisFlag,PathName1,FileNames1, PathName2, FileNames2, AnalysisParam);
    end
end


% --- Executes on button press in TrackDimButton.
function TrackDimButton_Callback(hObject, eventdata, handles)
% hObject    handle to TrackDimButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.PreAnalysis = [];
handles.Analysis = [];
handles.Parameters.Used.discardCheck = handles.Parameters.discardCheck;

if ~isfield(handles.Tracking, 'Tracks')
    errordlg('You need to track the particles first!');
    return
end;



if isfield(handles.Tracking, 'CheckTracks')
    
    button = questdlg(['Some tracks have already been edited. ';...
        'Do you want to continue or start over?'],'Edited Tracks Found', ...
        'Continue', 'Start Over','Continue');
    
else
    button = 'Start Over';
    
end

switch button
    
    case 'Start Over'
        TracksTemp = discardTracks(handles.Tracking.Tracks, ...
            handles.Parameters.discardCheck);
        if ~isfield(handles,'Process')
            handles.Process.ROIpos = [];
        else
            if ~isfield(handles.Process,'ROIpos')
                handles.Process.ROIpos = [];
            end
        end
        %         AnalysisType = questdlg('Do you want to adjust the individual Track Points, or just define the Track length & width?','Type of Checking','Points','Dimensions','Points');
        AnalysisType = 'Dimensions';
        if handles.isFitPSF
            if strcmp(AnalysisType,'Points')
                hCheckGUI = checkTracks2_GUI...
                    (handles.Data.imageStack, TracksTemp, handles.Tracking.Particles, handles.Process);
                %                 [handles.Tracking.CheckTracks handles.Tracking.CheckParticles] = checkTracks2_GUI...
                %                     (handles.Data.imageStack, TracksTemp, handles.Tracking.Particles, handles.Process);
                waitfor(hCheckGUI);
                handles.Tracking.CheckTracks = getappdata(0,'outTracks');
                handles.Tracking.CheckParticles = getappdata(0,'outParticles');
                handles.Tracking.TrackDim = [];
            else
                handles.Tracking.TrackDim = zeros(max(TracksTemp(:,4)),9);
                hCheckGUI = checkTracks3_GUI...
                    (handles.Data.imageStack, TracksTemp, handles.Tracking.Particles, handles.Process,handles.Tracking.TrackDim);
                %                 [handles.Tracking.CheckTracks handles.Tracking.CheckParticles,handles.Tracking.TrackDim] = checkTracks3_GUI...
                %                     (handles.Data.imageStack, TracksTemp, handles.Tracking.Particles, handles.Process,handles.Tracking.TrackDim);
                waitfor(hCheckGUI);
                handles.Tracking.CheckTracks = getappdata(0,'outTracks');
                handles.Tracking.CheckParticles = getappdata(0,'outParticles');
                handles.Tracking.TrackDim = getappdata(0,'outTrackDim');
                handles.Tracking.TrackDim(handles.Tracking.TrackDim(:,1) == 0,:) = [];
            end
        else
            if strcmp(AnalysisType,'Points')
                hCheckGUI = checkTracks2_GUI...
                    (handles.Data.imageStack, TracksTemp, handles.Tracking.Centroids, handles.Process);
                %                 [handles.Tracking.CheckTracks handles.Tracking.CheckParticles] = checkTracks2_GUI...
                %                     (handles.Data.imageStack, TracksTemp, handles.Tracking.Centroids, handles.Process);
                waitfor(hCheckGUI);
                handles.Tracking.CheckTracks = getappdata(0,'outTracks');
                handles.Tracking.CheckParticles = getappdata(0,'outParticles');
                handles.Tracking.TrackDim = [];
            else
                handles.Tracking.TrackDim = zeros(max(TracksTemp(:,4)),9);
                hCheckGUI = checkTracks3_GUI...
                    (handles.Data.imageStack, TracksTemp, handles.Tracking.Centroids, handles.Process,handles.Tracking.TrackDim);
                %                 [handles.Tracking.CheckTracks handles.Tracking.CheckParticles,handles.Tracking.TrackDim] = checkTracks3_GUI...
                %                     (handles.Data.imageStack, TracksTemp, handles.Tracking.Centroids, handles.Process,handles.Tracking.TrackDim);
                waitfor(hCheckGUI);
                handles.Tracking.CheckTracks = getappdata(0,'outTracks');
                handles.Tracking.CheckParticles = getappdata(0,'outParticles');
                handles.Tracking.TrackDim = getappdata(0,'outTrackDim');
                handles.Tracking.TrackDim(handles.Tracking.TrackDim(:,1) == 0,:) = [];
            end
        end
        
    case 'Continue'
        AnalysisType = 'Dimensions';
        if isfield(handles.Tracking, 'CheckParticles')
            if strcmp(AnalysisType,'Points')
                hCheckGUI = checkTracks2_GUI...
                    (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.CheckParticles, handles.Process);
                %                 [handles.Tracking.CheckTracks handles.Tracking.CheckParticles] = checkTracks2_GUI...
                %                     (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.CheckParticles, handles.Process);
                waitfor(hCheckGUI);
                handles.Tracking.CheckTracks = getappdata(0,'outTracks');
                handles.Tracking.CheckParticles = getappdata(0,'outParticles');
                handles.Tracking.TrackDim = [];
            else
                hCheckGUI = checkTracks3_GUI...
                    (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.CheckParticles, handles.Process,handles.Tracking.TrackDim);
                %                 [handles.Tracking.CheckTracks handles.Tracking.CheckParticles,handles.Tracking.TrackDim] = checkTracks3_GUI...
                %                     (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.CheckParticles, handles.Process,handles.Tracking.TrackDim);
                waitfor(hCheckGUI);
                handles.Tracking.CheckTracks = getappdata(0,'outTracks');
                handles.Tracking.CheckParticles = getappdata(0,'outParticles');
                handles.Tracking.TrackDim = getappdata(0,'outTrackDim');
                handles.Tracking.TrackDim(handles.Tracking.TrackDim(:,1) == 0,:) = [];
            end
        else
            
            if handles.isFitPSF
                if strcmp(AnalysisType,'Points')
                    hCheckGUI = checkTracks2_GUI...
                        (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.Particles, handles.Process);
                    %                     [handles.Tracking.CheckTracks handles.Tracking.CheckParticles] = checkTracks2_GUI...
                    %                         (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.Particles, handles.Process);
                    waitfor(hCheckGUI);
                    handles.Tracking.CheckTracks = getappdata(0,'outTracks');
                    handles.Tracking.CheckParticles = getappdata(0,'outParticles');
                    handles.Tracking.TrackDim = [];
                else
                    hCheckGUI = checkTracks3_GUI...
                        (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.Particles, handles.Process,handles.Tracking.TrackDim);
                    %                     [handles.Tracking.CheckTracks handles.Tracking.CheckParticles,handles.Tracking.TrackDim] = checkTracks3_GUI...
                    %                         (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.Particles, handles.Process,handles.Tracking.TrackDim);
                    waitfor(hCheckGUI);
                    handles.Tracking.CheckTracks = getappdata(0,'outTracks');
                    handles.Tracking.CheckParticles = getappdata(0,'outParticles');
                    handles.Tracking.TrackDim = getappdata(0,'outTrackDim');
                    handles.Tracking.TrackDim(handles.Tracking.TrackDim(:,1) == 0,:) = [];
                    
                end
            else
                if strcmp(AnalysisType,'Points')
                    hCheckGUI = checkTracks2_GUI...
                        (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.Centroids, handles.Process);
                    %                     [handles.Tracking.CheckTracks handles.Tracking.CheckParticles] = checkTracks2_GUI...
                    %                         (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.Centroids, handles.Process);
                    waitfor(hCheckGUI);
                    handles.Tracking.CheckTracks = getappdata(0,'outTracks');
                    handles.Tracking.CheckParticles = getappdata(0,'outParticles');
                    
                    handles.Tracking.TrackDim = [];
                else
                    hCheckGUI = checkTracks3_GUI...
                        (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.Centroids, handles.Process,handles.Tracking.TrackDim);
                    %                     [handles.Tracking.CheckTracks handles.Tracking.CheckParticles,handles.Tracking.TrackDim] = checkTracks3_GUI...
                    %                         (handles.Data.imageStack, handles.Tracking.CheckTracks, handles.Tracking.Centroids, handles.Process,handles.Tracking.TrackDim);
                    waitfor(hCheckGUI);
                    handles.Tracking.CheckTracks = getappdata(0,'outTracks');
                    handles.Tracking.CheckParticles = getappdata(0,'outParticles');
                    handles.Tracking.TrackDim = getappdata(0,'outTrackDim');
                    handles.Tracking.TrackDim(handles.Tracking.TrackDim(:,1) == 0,:) = [];
                end
            end
        end
end



% Enable the otpions to visualize hand checked tracks
set(handles.TrackPopUp, 'Value', 2);
set(handles.TrackPopUp, 'Enable', 'on');
handles.Current.Tracks = handles.Tracking.CheckTracks;
plotTracks(handles.Current.Tracks,handles.Current.Ix);


% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function Gmask_tracks_Callback(hObject, eventdata, handles)
% hObject    handle to Gmask_tracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

defPar(1) = 0.264;
defPar(2) = 1.46;
defPar(3) = 0.51;
defPar(4) = 7;

defaults = {num2str(defPar(1)), num2str(defPar(2)), ...
    num2str(defPar(3)), num2str(defPar(4))};
prompt = {'Pixel Size',...
    'Objective N.A.',...
    'Emission Wavelength',...
    'Window Size'};
dlgtitle = 'Parameters for Gaussian Fitting';
AnalysisParam = inputdlg(prompt,dlgtitle, 1, defaults);
if ~isempty(AnalysisParam)
    Param.pixelSize = str2double(AnalysisParam{1});
    Param.NA = str2double(AnalysisParam{2});
    Param.wavelength = str2double(AnalysisParam{3});
    Param.windowSize = str2double(AnalysisParam{4});
    
    if isfield(handles.Tracking,'CheckTracks')
        Tracks = handles.Tracking.CheckTracks;
    else
        Tracks = handles.Tracking.Tracks;
    end
    gmask_props = Particle_Intensities(handles.Data.imageStack,Tracks,Param);
    
    if isfield(handles,'PreAnalysis') && ~isempty(handles.PreAnalysis)
        handles.PreAnalysis.Tracks_um(:,8) = gmask_props(:,3);
    end
end

writeTrkFiles = questdlg('Do you want to create .trk files for each track?','Save .trk files','Yes','No','Yes');

if strcmp(writeTrkFiles,'Yes');
    baseName = [handles.Data.pathName filesep handles.Data.fileName(1:end-4)];
    for i = 1:max(handles.PreAnalysis.Tracks_um(:,4))
        d = num2str(i);
        curTrack = handles.Tracking.CheckTracks(handles.Tracking.CheckTracks(:,4) == i,:);
        curTrack2 = handles.PreAnalysis.Tracks_um(handles.PreAnalysis.Tracks_um(:,4) == i,:);
        fid = fopen([baseName, '_p',d, '.trk'],'w');
        for j = 1:size(curTrack,1)
            fprintf(fid,'%11.3f\t %11.3f\t %11.3f\t %11.3f\t %11.3f\n',curTrack(j,1),curTrack(j,2),curTrack2(j,8),curTrack(j,3)-1,0);
        end
        fclose(fid);
        
    end
end

guidata(hObject,handles);




% --- Executes on button press in AdjCntr.
function AdjCntr_Callback(hObject, eventdata, handles)
% hObject    handle to AdjCntr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.imagePopup,'Value') == 1
    handles.Current.clims = handles.Data.clims;
else
    handles.Current.clims = handles.Process.clims;
end
%update Intensity sliders & edit boxes
set(handles.BlackValSlider,'Val',handles.Current.clims(1));
set(handles.WhiteValSlider,'Val',handles.Current.clims(2));

set(handles.BlackValEdit,'String',num2str(round(double(handles.Current.clims(1))*100)/100));
set(handles.WhiteValEdit,'String',num2str(round(double(handles.Current.clims(2))*100)/100));

set(handles.axes1,'NextPlot','replacechildren');
imagesc(handles.Current.Image, handles.Current.clims);       % plot image
% axis image;
if get(handles.showParticles, 'Value')
    if handles.isFitPSF
        plotParticle(handles.Tracking.Particles, ...
            handles.Current.Ix, handles.isFitPSF);
    else
        plotParticle(handles.Tracking.Centroids, ...
            handles.Current.Ix, handles.isFitPSF);
    end
end

if get(handles.showTracks, 'Value')
    plotTracks(handles.Current.Tracks, handles.Current.Ix);
end

if isfield(handles.Process,'ROIpos')
    plotROI(handles.Process.ROIpos);
    
end
guidata(hObject,handles);

% --- Executes on slider movement.
function WhiteValSlider_Callback(hObject, eventdata, handles)
% hObject    handle to WhiteValSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
BlackLevel = get(handles.BlackValSlider,'Val');
WhiteLevel = get(hObject,'Val');
if WhiteLevel == 0
    WhiteLevel = 0.1;
end

if WhiteLevel <= BlackLevel
    BlackLevel = WhiteLevel-0.1;
    set(handles.BlackValSlider,'Val',BlackLevel);
end

set(handles.BlackValEdit,'String',num2str(round(double(BlackLevel)*100)/100));
set(handles.WhiteValEdit,'String',num2str(round(double(WhiteLevel)*100)/100));
handles.Current.clims = [BlackLevel WhiteLevel];

set(handles.axes1,'NextPlot','replacechildren');
imagesc(handles.Current.Image, handles.Current.clims);       % plot image
% axis image;
%Plot ROI

if get(handles.showParticles, 'Value')
    if handles.isFitPSF
        plotParticle(handles.Tracking.Particles, ...
            handles.Current.Ix, handles.isFitPSF);
    else
        plotParticle(handles.Tracking.Centroids, ...
            handles.Current.Ix, handles.isFitPSF);
    end
end

if get(handles.showTracks, 'Value')
    plotTracks(handles.Current.Tracks, handles.Current.Ix);
end

if isfield(handles.Process,'ROIpos')
    plotROI(handles.Process.ROIpos);
    
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function WhiteValSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WhiteValSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function BlackValSlider_Callback(hObject, eventdata, handles)
% hObject    handle to BlackValSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
MaxBrightnessStack = get(handles.BlackValSlider,'Max');

BlackLevel = get(hObject,'Val');
WhiteLevel = get(handles.WhiteValSlider,'Val');

if BlackLevel == MaxBrightnessStack
    BlackLevel = MaxBrightnessStack -0.1;
end
if BlackLevel >= WhiteLevel
    WhiteLevel = BlackLevel+0.1;
    set(handles.WhiteValSlider,'Val',WhiteLevel);
end

set(handles.BlackValEdit,'String',num2str(round(double(BlackLevel)*100)/100));
set(handles.WhiteValEdit,'String',num2str(round(double(WhiteLevel)*100)/100));
handles.Current.clims = [BlackLevel WhiteLevel];


set(handles.axes1,'NextPlot','replacechildren');
imagesc(handles.Current.Image, handles.Current.clims);       % plot image
% axis image;
if get(handles.showParticles, 'Value')
    if handles.isFitPSF
        plotParticle(handles.Tracking.Particles, ...
            handles.Current.Ix, handles.isFitPSF);
    else
        plotParticle(handles.Tracking.Centroids, ...
            handles.Current.Ix, handles.isFitPSF);
    end
end

if get(handles.showTracks, 'Value')
    plotTracks(handles.Current.Tracks, handles.Current.Ix);
end

if isfield(handles.Process,'ROIpos')
    plotROI(handles.Process.ROIpos);
    
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function BlackValSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BlackValSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function BlackValEdit_Callback(hObject, eventdata, handles)
% hObject    handle to BlackValEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BlackValEdit as text
%        str2double(get(hObject,'String')) returns contents of BlackValEdit as a double
MaxBrightnessStack = get(handles.BlackValSlider,'Max');

BlackLevel = str2double(get(hObject,'String'));
%Make sure we don't try to set the brightness outside of the slider limits
if BlackLevel < 0
    BlackLevel = 0;
elseif BlackLevel > MaxBrightnessStack
    BlackLevel = MaxBrightnessStack-0.1;
end

set(hObject,'String',num2str(round(double(BlackLevel)*100)/100));
WhiteLevel = str2double(get(handles.WhiteValEdit,'String'));

if BlackLevel >= WhiteLevel
    WhiteLevel = BlackLevel+0.1;
    set(handles.WhiteValEdit,'String',num2str(round(double(WhiteLevel)*100)/100));
end

set(handles.BlackValSlider,'Val',BlackLevel);
set(handles.WhiteValSlider,'Val',WhiteLevel);
handles.Current.clims = [BlackLevel WhiteLevel];


set(handles.axes1,'NextPlot','replacechildren');
imagesc(handles.Current.Image, handles.Current.clims);       % plot image
% axis image;
%Plot ROI

if get(handles.showParticles, 'Value')
    if handles.isFitPSF
        plotParticle(handles.Tracking.Particles, ...
            handles.Current.Ix, handles.isFitPSF);
    else
        plotParticle(handles.Tracking.Centroids, ...
            handles.Current.Ix, handles.isFitPSF);
    end
end

if get(handles.showTracks, 'Value')
    plotTracks(handles.Current.Tracks, handles.Current.Ix);
end

if isfield(handles.Process,'ROIpos')
    plotROI(handles.Process.ROIpos);
    
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function BlackValEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BlackValEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function WhiteValEdit_Callback(hObject, eventdata, handles)
% hObject    handle to WhiteValEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WhiteValEdit as text
%        str2double(get(hObject,'String')) returns contents of WhiteValEdit as a double
MaxBrightnessStack = get(handles.BlackValSlider,'Max');

WhiteLevel = str2double(get(hObject,'String'));
%Make sure we don't try to set the brightness outside of the slider limits
if WhiteLevel < 0
    WhiteLevel = 0.1;
elseif WhiteLevel > MaxBrightnessStack
    WhiteLevel = MaxBrightnessStack;
end

set(hObject,'String',num2str(round(double(WhiteLevel)*100)/100));
BlackLevel = str2double(get(handles.BlackValEdit,'String'));

if WhiteLevel <= BlackLevel
    BlackLevel = WhiteLevel-0.1;
    set(handles.BlackValEdit,'String',num2str(round(double(BlackLevel)*100)/100));
end

set(handles.BlackValSlider,'Val',BlackLevel);
set(handles.WhiteValSlider,'Val',WhiteLevel);
handles.Current.clims = [BlackLevel WhiteLevel];

set(handles.axes1,'NextPlot','replacechildren');
imagesc(handles.Current.Image, handles.Current.clims);       % plot image
% axis image;
if get(handles.showParticles, 'Value')
    if handles.isFitPSF
        plotParticle(handles.Tracking.Particles, ...
            handles.Current.Ix, handles.isFitPSF);
    else
        plotParticle(handles.Tracking.Centroids, ...
            handles.Current.Ix, handles.isFitPSF);
    end
end

if get(handles.showTracks, 'Value')
    plotTracks(handles.Current.Tracks, handles.Current.Ix);
end

if isfield(handles.Process,'ROIpos')
    plotROI(handles.Process.ROIpos);
    
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function WhiteValEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WhiteValEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function ImportMMorphROIs_Callback(hObject, eventdata, handles)
% hObject    handle to ImportMMorphROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[roiFile,filepath] = uigetfile('*.log','Open Metamorph ROI log file',handles.Data.pathName);

if ~iscell(roiFile)
    roiFile2{1,:} = roiFile;
    roiFile = roiFile2;
end

if roiFile{1,:} ~= 0
    imageStr = inputdlg('Enter the pattern to search for in the Outline file','File pattern',1,{handles.Data.fileName(1:7)});
    
    fid = fopen([filepath,filesep,roiFile{1,:}],'r');
    ind = 1;
    str_ind = 1;
    tmp = fscanf(fid,'%s',1);
    while ~strcmp(tmp,'')
        
        if strcmp(tmp,'Calibrations:')
            for m = 1:19
                tmp = fscanf(fid,'%s',1);
            end
        end
        
        if rem(ind,7) == 1
            Im_names{str_ind,:} = tmp;
            ind = ind + 1;
        elseif rem(ind,7) == 2
            RegName{str_ind,:} = tmp;
            ind = ind + 1;
        elseif rem(ind,7) == 3
            Area(str_ind,:) = str2double(tmp);
            ind = ind + 1;
        elseif rem(ind,7) == 4
            Left(str_ind,:) = str2double(tmp);
            ind = ind + 1;
        elseif rem(ind,7) == 5
            Top(str_ind,:) = str2double(tmp);
            ind = ind + 1;
        elseif rem(ind,7) == 6
            Width(str_ind,:) = str2double(tmp);
            ind = ind + 1;
        else
            Height(str_ind,:) = str2double(tmp);
            ind = ind + 1;
            str_ind = str_ind + 1;
        end
        tmp = fscanf(fid,'%s',1);
    end
    
    fclose(fid);
    
    %Determine the image names within the file
    im_name = unique(Im_names);
    for i = 1:length(im_name)
        str_ind_find = strfind(im_name,imageStr{1,:});
        if ~isempty(str_ind_find)
            str_ind2 = i;
            break
        end
    end
    im_ind = find(strcmp(im_name{str_ind2,:},Im_names));
    ImROIs.Left = Left(im_ind,:);
    ImROIs.Top = Top(im_ind,:);
    ImROIs.Width = Width(im_ind,:);
    ImROIs.Height = Height(im_ind,:);
    radii = ImROIs.Width/2;
    centers = [(ImROIs.Left + (ImROIs.Width/2)),(ImROIs.Top + (ImROIs.Width/2))];
    hold on
    ROI_list = cell(length(im_ind),1);
    for j = 1:length(im_ind)
        hCirc = viscircles(centers(j,:),radii(j,:));
        CircData = get(hCirc,'Children');
        x = CircData.XData;
        x = round(x);
        y = CircData.YData;
        y = round(y);
        y(isnan(x)) = [];
        x(isnan(x)) = [];
        x(isnan(y)) = [];
        y(isnan(y)) = [];
        handles.Process.ROIpos{j,:} = [x',y'];
        
        ROI_list{j,:} = ['ROI_',num2str(j)];
        [height,width] = size(handles.Data.imageStack(1).data);
        
        [meshX, meshY] = meshgrid(1:width, 1:height);
        ROIimage_tmp = zeros(width,height);
        CheckMatrix = (centers(j,1) - meshX).^2 + (centers(j,2)-meshY).^2;
        idx = find (CheckMatrix<= radii(j,:).^2);
        ROIimage_tmp(idx) = 1;
        handles.Process.ROIimage{j,:} = ROIimage_tmp;
        
        
    end
    hold off
    set(handles.ROIList,'String',ROI_list);
    handles.Process.ROIlabel = ROI_list;
    set(handles.axes1,'NextPlot','replacechildren');
    imagesc(handles.Current.Image,handles.Current.clims);
    %     axis image;
    plotROI(handles.Process.ROIpos);
    guidata(hObject,handles);
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SpltAlignMenu_Callback(hObject, eventdata, handles)
% hObject    handle to SpltAlignMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SplitAndAlignImages_GUI


% --------------------------------------------------------------------
function SaveFilterMenu_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFilterMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles.Process,'filterStack')
    if ~isempty(handles.Process.filterStack)
        
        [fname, pname] = uiputfile('*.tif','Save Filtered Image Stack',[handles.Data.pathName, handles.Data.fileName(1:end-4), '_filtered.tif']);
        
        if fname ~=0
            for i = 1:handles.Data.nImages
                imwrite(uint16(handles.Process.filterStack(i).data),fullfile(pname,fname),'WriteMode','append');
            end
            
        end
    end
end


% --------------------------------------------------------------------
function BatchRetrackMenu_Callback(hObject, eventdata, handles)
% hObject    handle to BatchRetrackMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
defPar(1) = str2double(get(handles.thresholdEdit,'String'));
defPar(2) = str2double(get(handles.windowEdit,'String'));

defPar(3) = str2double(get(handles.jumpEdit,'String'));
defPar(4) = str2double(get(handles.shTrackEdit,'String'));
defPar(5) = str2double(get(handles.gapsCheck,'String'));

defaults = {num2str(defPar(1)), num2str(defPar(2)), ...
    num2str(defPar(3)), num2str(defPar(4)), num2str(defPar(5))};
prompt = {'Threshold (0 = No Particle identification)',...
    'Window Size (0 = No Particle identification)',...
    'Maximum Jump (0 = Keep existing value)',...
    'Shortest Track(0 = Keep existing value)',...
    'Gaps to close(-1 = Keep existing value)'};
dlgtitle = 'Parameters for Batch Tracking';
AnalysisParam = inputdlg(prompt,dlgtitle, 1, defaults);

if isfield(handles,'Data')
    if isfield(handles.Data,'pathName')
        
        curdir = pwd;
        cd(handles.Data.pathName);
    end
end
hpass = str2double(get(handles.hpEdit,'String'));

ReTrack(hpass,str2double(AnalysisParam{1}),str2double(AnalysisParam{2}),str2double(AnalysisParam{3}),str2double(AnalysisParam{4}),str2double(AnalysisParam{5}));

if isfield(handles,'Data')
    if isfield(handles.Data,'pathName')
        
        cd(curdir);
    end
end



% --------------------------------------------------------------------
function NminCalcMenu_Callback(hObject, eventdata, handles)
% hObject    handle to NminCalcMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CalcNmin_GUI();


% --- Executes on button press in CopyROI.
function CopyROI_Callback(hObject, eventdata, handles)
% hObject    handle to CopyROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Determine state of buttons
showPart_select = get(handles.showParticles,'Value');
showPart_enable = get(handles.showParticles,'Enable');

% Disable the options to visualize localized particles
set(handles.showParticles, 'Value', 0);
set(handles.showParticles, 'Enable', 'off');
% delete(findobj(gca,'Color','r'));

% Disable the ROI buttons
set(handles.roiButton, 'Enable', 'off');
%get the current state of the other ROI buttons
ROIrem_enable = get(handles.roiRemovePush,'Enable');
ROIname_enable = get(handles.roiNamePush,'Enable');

set(handles.roiRemovePush,'Enable','off');
set(handles.roiNamePush,'Enable','off');

TrackPop_enable = get(handles.showTracks,'Enable');
% Disable the otpions to visualize hand checked tracks
% set(handles.TrackPopUp, 'Value',1);
set(handles.TrackPopUp, 'Enable', 'off');

showTrack_select = get(handles.showTracks,'Value');
showTrack_enable = get(handles.showTracks,'Enable');
% Disable the options to visualize the tracks
set(handles.showTracks, 'Value', 0);
set(handles.showTracks, 'Enable', 'off');
delete(findobj(gca,'Color','g'));

% Ask if you want to load a reference image for selecting the ROI
if ~isfield(handles.Process,'RefImage') || isempty(handles.Process.RefImage)
    answer = questdlg('Do you want to open a reference image');
else
    answer = 'Yes';
end

switch answer
    case 'Yes'
        %Open the reference image
        if ~isfield(handles.Process,'RefImage') || isempty(handles.Process.RefImage)
            [fileName, pathName]= ...
                uigetfile ('.tif','Open Reference Image', handles.Data.pathName);
            if fileName == 0
                set(handles.roiButton, 'Enable', 'on');
                return
            else
                [StackRef_stack, nImages] = TIFread([pathName, fileName]);
                RefDesc = 'Sum of Other';
                changeDesc_flag = 1;
            end
            
            
            
            
        else
            StackRef_stack(1).data = handles.Process.RefImage;
            nImages = 1;
            RefDesc = get(handles.RefImageDescription,'String');
            changeDesc_flag = 0;
        end
        
        
    case 'No'
        %Calculate average projection of the images
        if get(handles.imagePopup,'Value') == 1
            StackRef_stack = handles.Data.imageStack;
            nImages = handles.Data.nImages;
            RefDesc = 'Sum of Original';
            changeDesc_flag = 1;
%             SumImage = double(handles.Data.imageStack(1).data);
%             for imageIx =  2:handles.Data.nImages;
%                 SumImage = SumImage +  double(handles.Data.imageStack(imageIx).data);
%             end;
        else
            StackRef_stack = handles.Process.filterStack;
            nImages = handles.Data.nImages;
            RefDesc = 'Sum of Filtered';
            changeDesc_flag = 1;
%             SumImage = double(handles.Process.filterStack(1).data);
%             for imageIx =  2:handles.Data.nImages;
%                 SumImage = SumImage +  double(handles.Process.filterStack(imageIx).data);
%             end;
        end
        
        
    case 'Cancel'
        % Enable the ROI button
        
        set(handles.roiButton, 'Enable', 'on');
        return
        
end
for imgIx = 1:nImages
    StackRef(:,:,imgIx) = double(StackRef_stack(imgIx).data);
end
if ~isfield(handles.Process,'RefImage') || isempty(handles.Process.RefImage)
    answer = questdlg('What type of projection do you want to use for the Reference Image','Projection Type','Sum','Maximum','Sum');
else
    answer = 'Sum';
end

switch answer
    case 'Sum'
        SumImage = sum(StackRef,3);
        
    case 'Maximum'
        SumImage = max(StackRef,[],3);
        if changeDesc_flag == 1;
            RefDesc = ['Max', RefDesc(4:end)];
        end
    
end
   set(handles.RefImageDescription,'String',RefDesc);
        

%Plot average projection of the stack
% if strcmp(answer,'No')
%     if max(max(SumImage))/min(min(SumImage(SumImage > 0))) > 6
%         projlims = [min(min(SumImage(SumImage > 0))) max(max(SumImage))];
%     elseif max(max(SumImage))/min(min(SumImage(SumImage > 0))) > 4
%         projlims = [min(min(SumImage(SumImage > 0))) max(max(SumImage))];
%     else
%         projlims = [min(min(SumImage(SumImage > 0))) max(max(SumImage))];
%     end
% else
projlims = [min(min(SumImage(SumImage > 0))) max(max(SumImage))];
% end

if isfield(handles.Analysis,'HistD')||isfield(handles.Analysis,'BoundF')
    rmDataAns = questdlg('Analysis data already exists that will be removed if you specify a new ROI. Proceed?','Existing Data');
    switch rmDataAns
        case 'Yes'
            handles.Analysis = [];
        case 'No'
            return
    end
    
    
end


%Plot average projection of the stack
set(handles.axes1,'NextPlot','replacechildren');
imagesc(SumImage, projlims);

if isfield(handles.Process,'ROIpos')
    nROIs = length(handles.Process.ROIpos);
    plotROI(handles.Process.ROIpos);
    
else
    nROIs = 0;
end
set(handles.StatusText,'String','Select the center of the new ROI');
drawnow;

ROI_sel = get(handles.ROIList,'Value');

ROI_pos = handles.Process.ROIpos{ROI_sel};
ROI_pos_cent_old = mean(ROI_pos);
ROI_pos_rel = ROI_pos - repmat(mean(ROI_pos),size(ROI_pos,1),1);

hPoint = impoint(handles.axes1);

newPos = getPosition(hPoint);
delete(hPoint);
dist_Trans = newPos - ROI_pos_cent_old;
ROI_pos_new = ROI_pos_rel + repmat(newPos,size(ROI_pos_rel,1),1);
mask_ind = find(handles.Process.ROIimage{ROI_sel,:} == 1);
mask_ind_new = mask_ind + round(dist_Trans(2))+round(dist_Trans(1))*size(handles.Process.ROIimage{ROI_sel,:},1);

if min(mask_ind_new(:)) < 1 || max(mask_ind_new(:)) > numel(handles.Process.ROIimage{ROI_sel,:})
    errordlg('New ROI is too close to the image edge','Invalid ROI placement');
else
    ROI_label_def = {['ROI ' num2str(nROIs+1)]};
    ROI_label = inputdlg('Create a label for the new ROI:','ROI label',1,ROI_label_def);
    
    if ~isempty(ROI_label)
        ROI_list = get(handles.ROIList,'String');
        if isempty(ROI_list)
            ROI_list = cell(1);
            ROI_list{1,:} = ROI_label{1,:};
        else
            ROI_list{nROIs+1,:} = ROI_label{1,:};
        end
        set(handles.ROIList,'String',ROI_list);
        
        handles.Process.ROIpos{nROIs+1,:} = ROI_pos_new;
        
        handles.Process.ROIimage{nROIs+1,:} = false(size(handles.Process.ROIimage{ROI_sel,:}));
        
        
        handles.Process.ROIimage{nROIs+1,:}(mask_ind_new) = 1;
        
        handles.Process.ROIlabel = get(handles.ROIList,'String');
        
        
        %         axis image;
        plotROI(handles.Process.ROIpos);
        
        drawnow;
        %Add classes if desired
        if isempty(handles.Process.ROIClass)
            
            AddClassAns = questdlg('Do you want to separate the ROIs into 2 or more classifications?', 'Classify ROIs','Yes','No','Cancel','No');
            if strcmp(AddClassAns,'Cancel')
                delete(hPoly);
                set(handles.roiButton, 'Enable', 'on');
                set(handles.roiRemovePush,'Enable',ROIrem_enable);
                set(handles.roiNamePush,'Enable',ROIname_enable);
            elseif strcmp(AddClassAns,'No')
                handles.Process.ROIClass = {0};
            else
                ClassName = inputdlg('Specify the name of the new class:','New Classification',1,{'Nucleus'});
                handles.Process.ROIClass = cell(size(handles.Process.ROIpos));
                handles.Process.ROIClass{end,:} = ClassName(1);
                handles.Process.AllROIClasses = ClassName(1);
                set(handles.ROIClassCur_text,'String',ClassName{1});
                set(handles.ChangeROIClass,'Enable','on');
            end
        elseif iscell(handles.Process.ROIClass{1,:})
            ROIstring = handles.Process.AllROIClasses;
            ROIstring{end+1,1} = 'New...';
            ROIidx = ROIClassChooseDlg(ROIstring);
            if ROIidx < size(ROIstring,1)
                handles.Process.ROIClass{end+1,:} = ROIstring{ROIidx,:};
            else
                ClassName = inputdlg('Specify the name of the new class:','New Classification',1,{'Nucleus'});
                %Need to verify that it is really new
                isnew = 1;
                for i = 1:size(handles.Process.AllROIClasses,1)
                    if strcmpi(ClassName{1},handles.Process.AllROIClasses{i,:})
                        isnew = 0;
                        break;
                    end
                end
                if isnew == 1
                    handles.Process.AllROIClasses{end+1,:} = ClassName{1};
                    handles.Process.ROIClass{end+1,:} = ClassName{1};
                else
                    handles.Process.ROIClass{end+1,:} = handles.Process.AllROIClasses{i,:};
                end
                
                
            end
            set(handles.ROIClassCur_text,'String',handles.Process.ROIClass{end,:});
            set(handles.ChangeROIClass,'Enable','on');
        else
            set(handles.ROIClassCur_text,'String','Not defined');
            set(handles.ChangeROIClass,'Enable','on');
        end
        set(handles.axes1,'NextPlot','replacechildren');
        imagesc(handles.Current.Image, handles.Current.clims);
        plotROI(handles.Process.ROIpos);
        
        %Reassign Particle, Centroid & Tracking data
        if isfield(handles.Tracking,'Particles')
            handles.Tracking.Particles = InsideROIcheck2(handles.Tracking.Particles,handles.Process.ROIimage,0);
        end
        
        if isfield(handles.Tracking,'Centroids')
            handles.Tracking.Centroids = InsideROIcheck2(handles.Tracking.Centroids,handles.Process.ROIimage,0);
        end
        
        if isfield(handles.Tracking,'Tracks')
            handles.Tracking.Tracks = InsideROIcheck2(handles.Tracking.Tracks,handles.Process.ROIimage,0);
        end
        
        if isfield(handles.Tracking,'CheckTracks')
            handles.Tracking.CheckTracks = InsideROIcheck2(handles.Tracking.CheckTracks,handles.Process.ROIimage,0);
        end
        
        if isfield(handles.Tracking,'CheckParticles')
            handles.Tracking.CheckParticles = InsideROIcheck2(handles.Tracking.CheckParticles,handles.Process.ROIimage,0);
        end
        
        if isfield(handles,'PreAnalysis')
            if isfield(handles.PreAnalysis,'Tracks_um')
                pixelSize = handles.Parameters.Acquisition.pixelSize;
                
                fileName = [handles.Data.pathName,handles.Data.fileName];
                if isfield(handles.Tracking,'CheckTracks')
                    if size(handles.Tracking.CheckTracks,1) == size(handles.PreAnalysis.Tracks_um,1)
                        handles.PreAnalysis.Tracks_um(:,5) = handles.Tracking.CheckTracks(:,5);
                    else
                        [Tracks_um, NParticles, IntensityHist] = preProcess_noGUI(handles.Tracking.CheckTracks,handles.Data.imageStack,...
                            handles.Tracking.CheckParticles, pixelSize, handles.Data.nImages, fileName,handles.Process.ROIpos);
                        handles.PreAnalysis.Tracks_um = Tracks_um;
                        handles.PreAnalysis.NParticles = NParticles;
                        handles.PreAnalysis.IntensityHist = IntensityHist;
                    end
                else
                    if size(handles.Tracking.Tracks,1) == size(handles.PreAnalysis.Tracks_um,1)
                        handles.PreAnalysis.Tracks_um(:,5) = handles.Tracking.Tracks(:,5);
                    else
                        [Tracks_um, NParticles, IntensityHist] = preProcess_noGUI(handles.Tracking.Tracks,handles.Data.imageStack,...
                            handles.Tracking.Particles, pixelSize, handles.Data.nImages, fileName,handles.Process.ROIpos);
                        handles.PreAnalysis.Tracks_um = Tracks_um;
                        handles.PreAnalysis.NParticles = NParticles;
                        handles.PreAnalysis.IntensityHist = IntensityHist;
                    end
                end
            end
        end
    end
end
% Enable the ROI buttons
set(handles.roiButton, 'Enable', 'on');
set(handles.roiRemovePush,'Enable','on');
set(handles.roiNamePush,'Enable','on');
set(handles.CopyROI,'Enable','on');

set(handles.showParticles, 'Enable', showPart_enable);
set(handles.showParticles,'Value',showPart_select);

set(handles.showTracks, 'Enable', showTrack_enable);
set(handles.showTracks,'Value',showTrack_select);

set(handles.showTracks,'Enable',TrackPop_enable);
if showPart_select == 1
    if handles.isFitPSF
        plotParticle(handles.Tracking.Particles, handles.Current.Ix, ...
            handles.isFitPSF);
    else
        plotParticle(handles.Tracking.Centroids, handles.Current.Ix, ...
            handles.isFitPSF);
    end
end

if showTrack_select == 1;
    plotTracks(handles.Current.Tracks, handles.Current.Ix);
end

set(handles.StatusText,'String','');
drawnow;


% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function Photobleach_merge_Callback(hObject, eventdata, handles)
% hObject    handle to Photobleach_merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PhotoBleachCells;

% --------------------------------------------------------------------
function BleachCompare_Callback(hObject, eventdata, handles)
% hObject    handle to BleachCompare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
compareBleaching;


% --- Executes on button press in ChangeROIClass.
function ChangeROIClass_Callback(hObject, eventdata, handles)
% hObject    handle to ChangeROIClass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.Process.ROIClass) || (~isempty(handles.Process.ROIClass{1,:}) && handles.Process.ROIClass{1,:}{1}(1) == 0)
    ReallyAdd = questdlg('No classifications exist for current dataset, Do you want to add them?','Add Classifications','OK', 'Cancel','Cancel');
    
else
    ReallyAdd = 'OK';
end

if strcmp(ReallyAdd,'OK')
    
    
    
    curSel = get(handles.ROIList,'Value');
    if isempty(handles.Process.ROIClass) || (~isempty(handles.Process.ROIClass{1,:}) && handles.Process.ROIClass{1,:}{1}(1) == 0)
        handles.Process.ROIClass = cell(size(handles.Process.ROIimage));
        ClassName = inputdlg('Specify the name of the new class:','New Classification',1,{'Nucleus'});
        handles.Process.ROIClass{curSel,:} = ClassName(1);
        handles.Process.AllROIClasses = ClassName(1);
        set(handles.ROIClassCur_text,'String',ClassName{1});
    else
        ROIstring = handles.Process.AllROIClasses;
        ROIstring{end+1,1} = 'New...';
        ButPos = get(hObject,'Position');
        figPos = get(gcf,'Position');
        XYPos(1) = (figPos(1) + ButPos(1)*figPos(3));
        XYPos(2) = (figPos(2) + ButPos(2)*figPos(4));
        ROIidx = ROIClassChooseDlg(ROIstring,XYPos);
        if ROIidx < size(ROIstring,1)
            handles.Process.ROIClass{curSel,:} = ROIstring{ROIidx,:};
        else
            ClassName = inputdlg('Specify the name of the new class:','New Classification',1,{'Nucleus'});
            handles.Process.ROIClass{curSel,:} = ClassName{1};
            %Need to verify that it is actually new
            isnew = 1;
            for i = 1:size(handles.Process.AllROIClasses,1)
                if strcmpi(ClassName{1},handles.Process.AllROIClasses{i,:})
                    isnew = 0;
                    break;
                end
            end
            if isnew == 1
                handles.Process.AllROIClasses{end+1,:} = ClassName{1};
            else
                handles.Process.ROIClass{curSel,:} = handles.Process.AllROIClasses{i,:};
            end
        end
    end
    set(handles.ROIClassCur_text,'String',handles.Process.ROIClass{curSel,:});
    guidata(hObject,handles);
end



% --------------------------------------------------------------------
function ResAnalOpen_Callback(hObject, eventdata, handles)
% hObject    handle to ResAnalOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AnalyzeResTimes_GUI();


% --------------------------------------------------------------------
function DiffCoeffMerge_Callback(hObject, eventdata, handles)
% hObject    handle to DiffCoeffMerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles, 'Data')
    [FileNames,PathName,FilterIndex] = uigetfile('.mat',...
        'Select Mat files for the analysis of the residence time',...
        [handles.Data.pathName, handles.Data.fileName],'MultiSelect','on');
else
    [FileNames,PathName,FilterIndex] = uigetfile('.mat',...
        'Select Mat files for the analysis of the residence time',...
        'MultiSelect','on');
end
if ~iscell(FileNames)
    if FileNames == 0
        FileNames = [];
    else
        FileNames = {FileNames};
    end
    
end

if ~isempty(FileNames)
    defParam = {num2str(handles.Parameters.Acquisition.frameTime), '25'};
    prompt = {'Frame Interval','Number of Bins (0 = auto)'};
    dlgtitle = 'Enter Parameters';
    Params = inputdlg(prompt,dlgtitle, 1, defParam);
    frameTime = str2double(Params{1});
    nBins = str2double(Params{2});
    DiffCoeffMerge(FileNames,PathName,frameTime,nBins);
    
    
end


% --------------------------------------------------------------------
function ResBoxplots_Callback(hObject, eventdata, handles)
% hObject    handle to ResBoxplots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[process_files,pnames] = uigetfile('*.mat','Select the Residence time analysis files that you want to compare','','MultiSelect','on');

if ~iscell(process_files)
    process_files = {process_files};
end
if process_files{1} ~= 0
    ResLabels = cell(length(process_files),1);
    ResLabels_all = [];
    allResTimes = [];
    for i = 1:length(process_files)
        clear curResLabels;
        def_name = {''};
        prompt = ['Enter the label for the file: ', process_files{i}];
        dlgTitle = 'Enter Labels';
        ResLabel1 = inputdlg(prompt,dlgTitle,1,def_name);
        ResLabels{i} = ResLabel1{1};
        IN = load(fullfile(pnames,process_files{i}));
        resHist = IN.Results.Hist.Res_PB(:,1:2);
        resHist(:,2) = resHist(:,2)*(resHist(2,1) - resHist(1,1));
        nTracks = IN.Results.totalTracks;
        resHist(:,2) = round(nTracks*resHist(:,2));
        curResTimes = [];
        for j = 1:size(resHist,1)
            curPoint = ones(resHist(j,2),1);
            curPoint = curPoint.*resHist(j,1);
            curResTimes = [curResTimes;curPoint];
        end
        for j = 1:size(curResTimes,1)
            curResLabels{j,1} = ResLabels{i};
        end
        allResTimes = [allResTimes;curResTimes];
        ResLabels_all = [ResLabels_all;curResLabels];
    end
    figure; boxplot(allResTimes,ResLabels_all);
    ylabel('Residence Time (s)');
    
    
end


% --- Executes on button press in StandardROIbutton.
function StandardROIbutton_Callback(hObject, eventdata, handles)
% hObject    handle to StandardROIbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Determine state of buttons
showPart_select = get(handles.showParticles,'Value');
showPart_enable = get(handles.showParticles,'Enable');

% Disable the options to visualize localized particles
set(handles.showParticles, 'Value', 0);
set(handles.showParticles, 'Enable', 'off');
% delete(findobj(gca,'Color','r'));

% Disable the ROI buttons
set(handles.roiButton, 'Enable', 'off');
%get the current state of the other ROI buttons
ROIrem_enable = get(handles.roiRemovePush,'Enable');
ROIname_enable = get(handles.roiNamePush,'Enable');
ROIcopy_enable = get(handles.CopyROI,'Enable');


set(handles.roiRemovePush,'Enable','off');
set(handles.roiNamePush,'Enable','off');

TrackPop_enable = get(handles.showTracks,'Enable');
% Disable the otpions to visualize hand checked tracks
% set(handles.TrackPopUp, 'Value',1);
set(handles.TrackPopUp, 'Enable', 'off');

showTrack_select = get(handles.showTracks,'Value');
showTrack_enable = get(handles.showTracks,'Enable');
% Disable the options to visualize the tracks
set(handles.showTracks, 'Value', 0);
set(handles.showTracks, 'Enable', 'off');
delete(findobj(gca,'Color','g'));

% Ask if you want to load a reference image for selecting the ROI
if ~isfield(handles.Process,'RefImage') || isempty(handles.Process.RefImage)
    answer = questdlg('Do you want to open a reference image');
else
    answer = 'Yes';
end

switch answer
    case 'Yes'
        %Open the reference image
        if ~isfield(handles.Process,'RefImage') || isempty(handles.Process.RefImage)
            [fileName, pathName]= ...
                uigetfile ('.tif','Open Reference Image', handles.Data.pathName);
            if fileName == 0
                set(handles.roiButton, 'Enable', 'on');
                return
            else
                [StackRef_stack, nImages] = TIFread([pathName, fileName]);
                RefDesc = 'Sum of Other';
                changeDesc_flag = 1;
            end
            
            
            
            
        else
            StackRef_stack(1).data = handles.Process.RefImage;
            nImages = 1;
            RefDesc = get(handles.RefImageDescription,'String');
            changeDesc_flag = 0;
        end
        
        
    case 'No'
        %Calculate average projection of the images
        if get(handles.imagePopup,'Value') == 1
            StackRef_stack = handles.Data.imageStack;
            nImages = handles.Data.nImages;
            RefDesc = 'Sum of Original';
            changeDesc_flag = 1;
%             SumImage = double(handles.Data.imageStack(1).data);
%             for imageIx =  2:handles.Data.nImages;
%                 SumImage = SumImage +  double(handles.Data.imageStack(imageIx).data);
%             end;
        else
            StackRef_stack = handles.Process.filterStack;
            nImages = handles.Data.nImages;
            RefDesc = 'Sum of Filtered';
            changeDesc_flag = 1;
%             SumImage = double(handles.Process.filterStack(1).data);
%             for imageIx =  2:handles.Data.nImages;
%                 SumImage = SumImage +  double(handles.Process.filterStack(imageIx).data);
%             end;
        end
        
        
    case 'Cancel'
        % Enable the ROI button
        
        set(handles.roiButton, 'Enable', 'on');
        return
        
end
for imgIx = 1:nImages
    StackRef(:,:,imgIx) = double(StackRef_stack(imgIx).data);
end
if ~isfield(handles.Process,'RefImage') || isempty(handles.Process.RefImage)
    answer = questdlg('What type of projection do you want to use for the Reference Image','Projection Type','Sum','Maximum','Sum');
else
    answer = 'Sum';
end

switch answer
    case 'Sum'
        SumImage = sum(StackRef,3);
        
    case 'Maximum'
        SumImage = max(StackRef,[],3);
        if changeDesc_flag == 1;
            RefDesc = ['Max', RefDesc(4:end)];
        end
   
end
   set(handles.RefImageDescription,'String',RefDesc);
   

%Plot average projection of the stack
% if strcmp(answer,'No')
%     if max(max(SumImage))/min(min(SumImage(SumImage > 0))) > 6
%         projlims = [min(min(SumImage(SumImage > 0))) max(max(SumImage))];
%     elseif max(max(SumImage))/min(min(SumImage(SumImage > 0))) > 4
%         projlims = [min(min(SumImage(SumImage > 0))) max(max(SumImage))];
%     else
%         projlims = [min(min(SumImage(SumImage > 0))) max(max(SumImage))];
%     end
% else
projlims = [min(min(SumImage(SumImage > 0))) max(max(SumImage))];
% end
set(handles.axes1,'NextPlot','replacechildren');
imagesc(SumImage, projlims);


%determine if the ROI type and size have been set
if isfield(handles.Parameters,'StandardROI')
    ROItype = handles.Parameters.StandardROI.ROItype;
    ROIsize = handles.Parameters.StandardROI.ROIsize;
    
    hSetROIParam = SetROIParamGUI(ROItype,ROIsize);
else
    hSetROIParam = SetROIParamGUI();
end

waitfor(hSetROIParam);
handles.Parameters.StandardROI.ROItype = getappdata(0,'ROItype');
handles.Parameters.StandardROI.ROIsize = getappdata(0,'ROIsize');



if isfield(handles.Analysis,'HistD')||isfield(handles.Analysis,'BoundF')
    rmDataAns = questdlg('Analysis data already exists that will be removed if you specify a new ROI. Proceed?','Existing Data');
    switch rmDataAns
        case 'Yes'
            handles.Analysis = [];
        case 'No'
            return
    end
    
    
end




if isfield(handles.Process,'ROIpos')
    nROIs = length(handles.Process.ROIpos);
    plotROI(handles.Process.ROIpos);
    
else
    nROIs = 0;
end
set(handles.StatusText,'String','Select the center of the new ROI');
drawnow;

hPoint = impoint(handles.axes1);

newPos = getPosition(hPoint);
delete(hPoint);

%Generate a default ROI centered at the selected point
newROIim = false(size(handles.Data.imageStack(1).data));

if handles.Parameters.StandardROI.ROItype == 1
    x1 = -handles.Parameters.StandardROI.ROIsize/2:0.1:handles.Parameters.StandardROI.ROIsize/2;
    x2 = -1*x1;
    
    y1 = sqrt((handles.Parameters.StandardROI.ROIsize/2)^2 - x1.^2);
    y2 = -1*y1;
    x = [x1'; x2'];
    
    y = [y1';y2'];
    x(end) = [];
    y(end) = [];
    x = x + newPos(1);
    y = y + newPos(2);
    
    
    
else
    x = [-handles.Parameters.StandardROI.ROIsize/2;  handles.Parameters.StandardROI.ROIsize/2; ...
        handles.Parameters.StandardROI.ROIsize/2; -handles.Parameters.StandardROI.ROIsize/2];
    y = [-handles.Parameters.StandardROI.ROIsize/2;  -handles.Parameters.StandardROI.ROIsize/2; ...
        handles.Parameters.StandardROI.ROIsize/2; handles.Parameters.StandardROI.ROIsize/2];
    
    x = x + newPos(1);
    y = y + newPos(2);
    
    
    
end
if min(x) < 1 || min(y) < 1 || max(x) > handles.Data.imageStack(1).width || max(y) > handles.Data.imageStack(1).height
    errordlg('New ROI is too close to the image edge','Invalid ROI placement');
    set(handles.roiRemovePush,'Enable',ROIrem_enable);
    set(handles.roiNamePush,'Enable',ROIname_enable);
    set(handles.CopyROI,'Enable',ROIcopy_enable);
else
    if handles.Parameters.StandardROI.ROItype == 1
        for i = round(min(x)):round(max(x))
            for j = round(min(y)):round(max(y));
                if sqrt((i - newPos(1)).^2 + (j - newPos(2)).^2) <= handles.Parameters.StandardROI.ROIsize/2
                    newROIim(j,i) = 1;
                end
            end
        end
    else
        for i = round(min(x)):round(max(x))
            for j = round(min(y)):round(max(y));
                newROIim(j,i) = 1;
            end
        end
    end
    ROI_label_def = {['ROI ' num2str(nROIs+1)]};
    ROI_label = inputdlg('Create a label for the new ROI:','ROI label',1,ROI_label_def);
    
    if ~isempty(ROI_label)
        ROI_list = get(handles.ROIList,'String');
        if isempty(ROI_list)
            ROI_list = cell(1);
            ROI_list{1,:} = ROI_label{1,:};
        else
            ROI_list{nROIs+1,:} = ROI_label{1,:};
        end
        set(handles.ROIList,'String',ROI_list);
        set(handles.ROIList,'Value',nROIs+1);
        
        handles.Process.ROIpos{nROIs+1,:} = [x,y];
        
        
        handles.Process.ROIimage{nROIs+1,:} = newROIim;
        handles.Process.ROIlabel = get(handles.ROIList,'String');
        plotROI(handles.Process.ROIpos);
        drawnow;
        handles.Process.RefImage = SumImage;
        set(handles.RemRefImage,'Enable','on');
        %Add classes if desired
        if isempty(handles.Process.ROIClass)
            
            AddClassAns = questdlg('Do you want to separate the ROIs into 2 or more classifications?', 'Classify ROIs','Yes','No','Cancel','No');
            if strcmp(AddClassAns,'Cancel')
                
                set(handles.roiButton, 'Enable', 'on');
                set(handles.roiRemovePush,'Enable',ROIrem_enable);
                set(handles.roiNamePush,'Enable',ROIname_enable);
            elseif strcmp(AddClassAns,'No')
                handles.Process.ROIClass = {0};
            else
                ClassName = inputdlg('Specify the name of the new class:','New Classification',1,{'Nucleus'});
                handles.Process.ROIClass = cell(size(handles.Process.ROIpos));
                handles.Process.ROIClass{end,:} = ClassName(1);
                handles.Process.AllROIClasses = ClassName(1);
                set(handles.ROIClassCur_text,'String',ClassName{1});
                set(handles.ChangeROIClass,'Enable','on');
            end
        elseif iscell(handles.Process.ROIClass{1,:})
            ROIstring = handles.Process.AllROIClasses;
            ROIstring{end+1,1} = 'New...';
            ROIidx = ROIClassChooseDlg(ROIstring);
            if ROIidx < size(ROIstring,1)
                handles.Process.ROIClass{end+1,:} = ROIstring{ROIidx,:};
            else
                ClassName = inputdlg('Specify the name of the new class:','New Classification',1,{'Nucleus'});
                %Need to verify that it is really new
                isnew = 1;
                for i = 1:size(handles.Process.AllROIClasses,1)
                    if strcmpi(ClassName{1},handles.Process.AllROIClasses{i,:})
                        isnew = 0;
                        break;
                    end
                end
                if isnew == 1
                    handles.Process.AllROIClasses{end+1,:} = ClassName{1};
                    handles.Process.ROIClass{end+1,:} = ClassName{1};
                else
                    handles.Process.ROIClass{end+1,:} = handles.Process.AllROIClasses{i,:};
                end
                
                
            end
            set(handles.ROIClassCur_text,'String',handles.Process.ROIClass{end,:});
            set(handles.ChangeROIClass,'Enable','on');
        else
            set(handles.ROIClassCur_text,'String','Not defined');
            set(handles.ChangeROIClass,'Enable','on');
        end
        set(handles.axes1,'NextPlot','replacechildren');
        imagesc(handles.Current.Image, handles.Current.clims);
        %     axis image;
        plotROI(handles.Process.ROIpos);
        set(handles.CopyROI,'Enable','on');
        set(handles.roiRemovePush,'Enable','on');
        set(handles.roiNamePush,'Enable','on');
        set(handles.CopyROI,'Enable','on');
        
    else
        set(handles.roiRemovePush,'Enable',ROIrem_enable);
        set(handles.roiNamePush,'Enable',ROIname_enable);
        set(handles.CopyROI,'Enable',ROIcopy_enable);
    end
end
set(handles.StatusText,'String','');
drawnow;

set(handles.roiButton, 'Enable', 'on');



guidata(hObject,handles);


% --- Executes on button press in RemRefImage.
function RemRefImage_Callback(hObject, eventdata, handles)
% hObject    handle to RemRefImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Process = rmfield(handles.Process,'RefImage');
set(hObject,'Enable','off');
set(handles.RefImageDescription,'String','N/A');

guidata(hObject,handles);


% --------------------------------------------------------------------
function SaveMatNoImage_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMatNoImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Assign .mat file where variables will be saved.
defName = [handles.Data.pathName, handles.Data.fileName];
SearchStr = '(.*)\.\w*';
defName = regexprep(defName, SearchStr, '$1');
FilterSpec = {'*.mat'};
[FileNameOut,PathNameOut] = uiputfile(FilterSpec,'Save Variables',defName);
if FileNameOut ~=0
    set(handles.StatusText,'String','Saving MAT file...');
    drawnow;
    Version = 2.0;
    
    %     save([PathNameOut, FileNameOut], 'Version');
    Results.Data.fileName = handles.Data.fileName;
    Results.Data.pathName = handles.Data.pathName;
    Results.Data.nImages = handles.Data.nImages;
    Results.Data.clims = handles.Data.clims;
    Results.Data.MaxBrightnessStack = handles.Data.MaxBrightnessStack;
    Results.Parameters = handles.Parameters.Used;
    Results.Process.ROIClass = handles.Process.ROIClass;
    Results.Process.clims = handles.Process.clims;
    Results.Process.MaxBrightnessStack = handles.Process.MaxBrightnessStack;
    Results.Process.ROIpos = handles.Process.ROIpos;
    Results.Process.ROIimage = handles.Process.ROIimage;
    Results.Process.ROIlabel = handles.Process.ROIlabel;
    Results.Process.RefImage = handles.Process.RefImage;
    Results.Tracking = handles.Tracking;
    Results.Analysis = handles.Analysis;
    Results.PreAnalysis = handles.PreAnalysis;
    Results.isFitPSF = handles.isFitPSF;
    save([PathNameOut, FileNameOut], 'Version','Results', '-v7.3');
    set(handles.StatusText,'String','Saving MAT file...Done');
    drawnow;
end
