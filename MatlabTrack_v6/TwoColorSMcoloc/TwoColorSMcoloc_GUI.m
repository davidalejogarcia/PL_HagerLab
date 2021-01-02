function varargout = TwoColorSMcoloc_GUI(varargin)
% TWOCOLORSMCOLOC_GUI MATLAB code for TwoColorSMcoloc_GUI.fig
%      TWOCOLORSMCOLOC_GUI, by itself, creates a new TWOCOLORSMCOLOC_GUI or raises the existing
%      singleton*.
%
%      H = TWOCOLORSMCOLOC_GUI returns the handle to a new TWOCOLORSMCOLOC_GUI or the handle to
%      the existing singleton*.
%
%      TWOCOLORSMCOLOC_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TWOCOLORSMCOLOC_GUI.M with the given input arguments.
%
%      TWOCOLORSMCOLOC_GUI('Property','Value',...) creates a new TWOCOLORSMCOLOC_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TwoColorSMcoloc_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TwoColorSMcoloc_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TwoColorSMcoloc_GUI

% Last Modified by GUIDE v2.5 28-Feb-2017 11:59:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @TwoColorSMcoloc_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @TwoColorSMcoloc_GUI_OutputFcn, ...
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


% --- Executes just before TwoColorSMcoloc_GUI is made visible.
function TwoColorSMcoloc_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TwoColorSMcoloc_GUI (see VARARGIN)

% Choose default command line output for TwoColorSMcoloc_GUI
handles.output = hObject;
DefaultBoundParam = load('IntegratedGUI_defaults.mat', ...
    'ThreshL','ThreshH', 'minBoundFrames');
handles.Parameters.Bound.ThreshL = DefaultBoundParam.ThreshL;
handles.Parameters.Bound.ThreshH = DefaultBoundParam.ThreshH;
handles.Parameters.Bound.minBoundFrames = DefaultBoundParam.minBoundFrames;

handles.isAlign = 0;
handles.PreAnalysis1 = [];
%Load in data if it was requested in Parent GUI
if ~isempty(varargin)
    %Set data structures
    
    handles.Data1 = varargin{1};
    handles.Process1 = varargin{2};
    handles.Tracking1= varargin{3};
    handles.PreAnalysis1 = varargin{4};
    if ~isempty(varargin{5})
        if isfield(varargin{5},'Acquisition')
            handles.Parameters.Used1.Acquisition = varargin{5}.Acquisition;
            handles.Parameters.Acquisition1 = varargin{5}.Acquisition;
            if isfield(varargin{5}.Acquisition,'frameTime')
                set(handles.FrameIntervalText,'String',num2str(varargin{5}.Acquisition.frameTime,3));
            else
                warndlg('It appears that this dataset was not preprocessed.', 'No Preprocess Data');
                
            end
            if isfield(varargin{5}.Acquisition,'pixelSize')
                set(handles.PxSizeText,'String',num2str(varargin{5}.Acquisition.pixelSize,3));
            end
        else
            warndlg('It appears that this dataset was not preprocessed.', 'No Preprocess Data');
            handles.Parameters.Acquisition1 = [];
        end
        if isfield(varargin{5},'discardCheck')
            handles.Parameters.Used1.discardCheck = varargin{5}.discardCheck;
        end
        if isfield(varargin{5},'Analysis');
            handles.Parameters.Used1.Analysis = varargin{5}.Analysis;
        end
        if isfield(varargin{5},'Bound');
            handles.Parameters.Used1.Bound = varargin{5}.Bound;
        end
        if isfield(varargin{5},'Used');
            if isfield(varargin{5}.Used,'Tracking')
                
                handles.Parameters.Used1.Tracking = varargin{5}.Used.Tracking;
            end
        end
    end
    
    % Enable some controls
    set(handles.showParticlesCheck,'Enable','on');
    set(handles.showTracksCheck,'Enable','on');
    set(handles.showTracksCheck,'Value',1);
    set(handles.LoadData2Push,'Enable','on');
    set(handles.imageSlider,'Enable','on','Min',1,'Max',handles.Data1.nImages,'Value',1,'SliderStep',[1/(handles.Data1.nImages-1) 1/(handles.Data1.nImages-1)]);
    
    
    %Disply the first image
    handles.Current.Ind = 1;
    set(handles.ImageCount,'String',['Image ' num2str(handles.Current.Ind) '/' num2str(handles.Data1.nImages)]);
    axes(handles.Channel1data);
    colormap gray;
    handles.Current.Image1 = ...
        handles.Data1.imageStack(handles.Current.Ind).data;
    %calculate the minimum & maximum values in the stack
    Mn_o = 100000;
    Mx_o = 0;
    
    for i = 1:handles.Data1.nImages
        if min(handles.Data1.imageStack(i).data(:)) < Mn_o
            Mn_o = min(handles.Data1.imageStack(i).data(:));
        end
        if max(handles.Data1.imageStack(i).data(:)) > Mx_o
            Mx_o = max(handles.Data1.imageStack(i).data(:));
        end
    end
    %     if Mx_o/Mn_o > 6
    %         handles.Current.clims1 = [Mn_o 0.5*Mx_o];
    %     elseif Mx_o/Mn_o > 4
    %         handles.Current.clims1 = [Mn_o 0.5*Mx_o];
    %     else
    handles.Current.clims1 = [Mn_o Mx_o];
    %     end
    handles.Data1.clims_orig = handles.Current.clims1;
    
    imagesc(handles.Current.Image1, handles.Current.clims1);
    axis image;
    
    %Enable other fields & plot if the data is available
    if isfield(handles.Process1,'filterStack')
        set(handles.ImageViewPop,'Enable','on');
        Mn_f = 100000;
        Mx_f = 0;
        for i = 1:handles.Data1.nImages
            if min(handles.Process1.filterStack(i).data(:)) < Mn_f
                Mn_f = min(handles.Process1.filterStack(i).data(:));
            end
            if max(handles.Process1.filterStack(i).data(:)) > Mx_f
                Mx_f = max(handles.Process1.filterStack(i).data(:));
            end
            %             if Mx_f/Mn_f > 6
            %                 handles.Data1.clims_filt = [Mn_f 0.25*Mx_f];
            %             elseif Mx_o/Mn_o > 4
            %                 handles.Data1.clims_filt = [Mn_f 0.5*Mx_f];
            %             else
            handles.Data1.clims_filt = [Mn_f Mx_f];
            %             end
            
        end
        
    end
    
    
    if isfield(handles.Process1,'ROIpos')
        set(handles.showROIsCheck,'Enable','on');
        set(handles.showROIsCheck,'Value',1);
        plotROI(handles.Process1.ROIpos);
    end
    
    if isfield(handles.Tracking1,'Particles')
        handles.isFitPSF1 = 1;
    else
        handles.isFitPSF1 = 0;
    end
    
    set(handles.showParticlesCheck,'Value',1)
    if handles.isFitPSF1 == 1
        plotParticle(handles.Tracking1.Particles,handles.Current.Ind,handles.isFitPSF1);
    else
        plotParticle(handles.Tracking1.Centroids,handles.Current.Ind,handles.isFitPSF1);
    end
    if isfield(handles.Tracking1,'CheckTracks')
        set(handles.TrackViewPop,'Enable','on');
        set(handles.TrackViewPop,'Value',2);
        handles.Current.Tracks1 = handles.Tracking1.CheckTracks;
    else
        set(handles.TrackViewPop,'Value',1);
        handles.Current.Tracks1 = handles.Tracking1.Tracks;
    end
    plotTracks(handles.Current.Tracks1, handles.Current.Ind);
    
    
    
    
    
    
    
end


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TwoColorSMcoloc_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TwoColorSMcoloc_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in showParticlesCheck.
function showParticlesCheck_Callback(hObject, eventdata, handles)
% hObject    handle to showParticlesCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showParticlesCheck

toShow = get(hObject,'Value');

if toShow == 1
    axes(handles.Channel1data);
    if handles.isFitPSF1
        plotParticle(handles.Tracking1.Particles,handles.Current.Ind,handles.isFitPSF1);
    else
        plotParticle(handles.Tracking1.Centroids,handles.Current.Ind,handles.isFitPSF1);
    end
    if isfield(handles,'Tracking2')
        axes(handles.Channel2data);
        if handles.isFitPSF2
            plotParticle(handles.Tracking2.Particles,handles.Current.Ind,handles.isFitPSF2);
        else
            plotParticle(handles.Tracking2.Centroids,handles.Current.Ind,handles.isFitPSF2);
        end
    end
    
else
    delete(findobj(handles.Channel1data,'MarkerSize',12));
    delete(findobj(handles.Channel2data,'MarkerSize',12));
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in showTracksCheck.
function showTracksCheck_Callback(hObject, eventdata, handles)
% hObject    handle to showTracksCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showTracksCheck
toShow = get(hObject,'Value');
if toShow == 1
    axes(handles.Channel1data);
    plotTracks(handles.Current.Tracks1, handles.Current.Ind);
    if isfield(handles.Current,'Tracks2')
        axes(handles.Channel2data);
        plotTracks(handles.Current.Tracks2, handles.Current.Ind);
    end
else
    delete(findobj(handles.Channel1data,'Color','g'));
    delete(findobj(handles.Channel2data,'Color','g'));
end;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in showROIsCheck.
function showROIsCheck_Callback(hObject, eventdata, handles)
% hObject    handle to showROIsCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showROIsCheck
toShow = get(hObject,'Value');

if toShow == 1
    axes(handles.Channel1data);
    plotROI(handles.Process1.ROIpos);
    if isfield(handles,'Process2')
        if isfield(handles.Process2,'ROIpos')
            axes(handles.Channel2data);
            plotROI(handles.Process2.ROIpos);
        end
    end
else
    delete(findobj(handles.Channel1data,'Linewidth',1));
    delete(findobj(handles.Channel2data,'Linewidth',1));
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in LoadData1Push.
function LoadData1Push_Callback(hObject, eventdata, handles)
% hObject    handle to LoadData1Push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName,FilterIndex] = ...
    uigetfile('.mat','Select Matlab File with tracked data');

if FileName ~=0;
    ListVar = whos('-file',[PathName,FileName]);
else
    return
end
% Reset the data content of the GUI
handles.Process1 = [];
handles.Parameters.Used1 = [];
handles.Tracking1 = [];
handles.Analysis1 = [];
handles.PreAnalysis1 = [];
handles.Current = [];
handles.Current.Ind = 1;
handles.isFitPSF1 = 0;

IN = load([PathName,FileName],'Version');

if IN.Version == 2.0;
    
    %Load all the contents of the structure.
    IN = load([PathName,FileName], 'Results');
    
    handles.Data1 =  IN.Results.Data;
%     handles.Data1.fileName = IN.Results.Data.fileName;
%     handles.Data1.pathName = IN.Results.Data.pathName;
    Parameters = IN.Results.Parameters;
    handles.Process1 = IN.Results.Process;
    handles.Tracking1 = IN.Results.Tracking;
    handles.isFitPSF1 = IN.Results.isFitPSF;
    handles.PreAnalysis1 = IN.Results.PreAnalysis;
    
    % Disable some of the controls
    % [They will be re-enabled if the variables are found within the input file]
    
    % Disable the otpions to visualize hand checked tracks
    set(handles.TrackViewPop, 'Value',1);
    set(handles.TrackViewPop, 'Enable', 'off');
    
    % Disable the options to visualize filtered images
    set(handles.ImageViewPop, 'Value',1);
    set(handles.ImageViewPop, 'Enable', 'off');
    
    % Disable the options to visualize localized particles
    set(handles.showParticlesCheck, 'Value', 0);
    set(handles.showParticlesCheck, 'Enable', 'off');
    
    % Disable the options to visualize the tracks
    set(handles.showTracksCheck, 'Value', 0);
    set(handles.showTracksCheck, 'Enable', 'off');
    
    %Enable filter button & roiButton
    %     set(handles.filterButton,'Enable','on');
    %     set(handles.roiButton,'Enable','on');
    %     set(handles.lpEdit,'Enable','on');
    %     set(handles.hpEdit,'Enable','on');
    
    
    % reset current index.
    handles.Current = [];
    handles.Current.Ind = 1;
    
    % Select current image and set colormap limits.
    handles.Current.Image1 = handles.Data1.imageStack(handles.Current.Ind).data;
    %calculate the minimum & maximum values in the stack
    Mn_o = 100000;
    Mx_o = 0;
    
    for i = 1:handles.Data1.nImages
        if min(handles.Data1.imageStack(i).data(:)) < Mn_o
            Mn_o = min(handles.Data1.imageStack(i).data(:));
        end
        if max(handles.Data1.imageStack(i).data(:)) > Mx_o
            Mx_o = max(handles.Data1.imageStack(i).data(:));
        end
    end
    
    handles.Current.clims1 = [Mn_o Mx_o];
    
    handles.Data1.clims_orig = handles.Current.clims1;
    set(handles.BlackSlider_C1,'Min',0,'Max',2*double(Mx_o),'SliderStep',[1/((2*double(Mx_o))-1),10/((2*double(Mx_o))-1)],'Value',double(Mn_o));
    set(handles.WhiteSlider_C1,'Min',0,'Max',2*double(Mx_o),'SliderStep',[1/((2*double(Mx_o))-1),10/((2*double(Mx_o))-1)],'Value',double(Mx_o));
    set(handles.BlackEdit_C1,'String',num2str(Mn_o));
    set(handles.WhiteEdit_C1,'String',num2str(Mx_o));
    
    %Enable other fields & plot if the data is available
    if isfield(handles.Process1,'filterStack')
        set(handles.ImageViewPop,'Enable','on');
        Mn_f = 100000;
        Mx_f = 0;
        for i = 1:handles.Data1.nImages
            if min(handles.Process1.filterStack(i).data(:)) < Mn_f
                Mn_f = min(handles.Process1.filterStack(i).data(:));
            end
            if max(handles.Process1.filterStack(i).data(:)) > Mx_f
                Mx_f = max(handles.Process1.filterStack(i).data(:));
            end
            %             if Mx_f/Mn_f > 6
            %                 handles.Data1.clims_filt = [Mn_f 0.25*Mx_f];
            %             elseif Mx_o/Mn_o > 4
            %                 handles.Data1.clims_filt = [Mn_f 0.5*Mx_f];
            %             else
            handles.Data1.clims_filt = [Mn_f Mx_f];
            %             end
            
        end
        
    end
    
    % Display current Image
    
    colormap(gray);
    axes(handles.Channel1data);
    imagesc(handles.Current.Image1, handles.Current.clims1);
    axis image;
    
    % Enable the image slider
    if handles.Data1.nImages > 1
        set(handles.imageSlider,'Enable','on');
        % Set Minimum of image slider
        set(handles.imageSlider,'Min', 1);
        % Set Maximum of image slider
        set(handles.imageSlider,'Max', handles.Data1.nImages);
        % Set Value of image slider
        set(handles.imageSlider,'Value', handles.Current.Ind);
        % Set Step of image slider
        set(handles.imageSlider, 'SliderStep', [1/(handles.Data1.nImages-1)...
            1/(handles.Data1.nImages-1)]);
        % Set Image counter
        set(handles.ImageCount, 'String', ['Image ', ...
            num2str(handles.Current.Ind),'/',num2str(handles.Data1.nImages)]);
    end;
    
    
    %     % Check if there is the filtered stack
    %     if isfield(handles.Process1, 'filterStack');
    %         % Enable the options to visualize filtered images
    %         set(handles.ImageViewPop, 'Value',1);
    %         set(handles.ImageViewPop, 'Enable', 'on');
    % %         set(handles.thresholdEdit,'Enable','on');
    % %         set(handles.windowEdit,'Enable','on');
    % %         set(handles.fitPSF,'Enable','on');
    % %         set(handles.findpeakButton,'Enable','on');
    %
    %
    %     end
    %Check if there is ROI data
    if isfield(handles.Process1,'ROIpos')
        nROIs1 = length(handles.Process1.ROIpos);
        set(handles.showROIsCheck,'Enable','on');
        set(handles.showROIsCheck,'Value',1);
        axes(handles.Channel1data);
        plotROI(handles.Process1.ROIpos);
        
        %         set(handles.roiRemovePush,'Enable','on');
        %         set(handles.roiNamePush,'Enable','on');
    end
    
    
    
    % Check if there is the data about the particles;
    if isfield(handles.Tracking1, 'Centroids')||isfield(handles.Tracking1,'Particles');
        % Check if there particles are fit with gaussians
        
        
        if isfield(handles.Tracking1, 'Particles');
            
            Centroids =  handles.Tracking1.Particles;
            
            
        else
            
            Centroids =  handles.Tracking1.Centroids;
            
        end
        axes(handles.Channel1data);
        plotParticle(Centroids,handles.Current.Ind,handles.isFitPSF1);
        
        set(handles.showParticlesCheck,'Enable','on');
        set(handles.showParticlesCheck,'Value',1);
        %         handles.arePeaksFiltered = 0;
        %         set(handles.Intensity_thr,'Enable','on');
        %         set(handles.sigmaLow,'Enable','on');
        %         set(handles.sigmaHigh,'Enable','on');
        %         set(handles.filterPeaks,'Enable','on');
        %         set(handles.editPeakButton,'Enable','on');
        %         set(handles.xPlot,'Enable','on');
        %         set(handles.yPlot,'Enable','on');
        %         set(handles.PlotVariables,'Enable','on');
        %         set(handles.jumpEdit,'Enable','on');
        %         set(handles.sigmaLow,'Enable','on');
        %         set(handles.shTrackEdit,'Enable','on');
        %         set(handles.gapsCheck,'Enable','on');
        %         set(handles.trackButton,'Enable','on');
        
    end
    
    % Check if there is information about tracks;
    if isfield(handles.Tracking1, 'Tracks')
        % Check if there are hand checked tracks
        if isfield(handles.Tracking1, 'CheckTracks')
            handles.Current.Tracks1 = handles.Tracking1.CheckTracks;
            set(handles.TrackViewPop, 'Value',2);
            set(handles.TrackViewPop, 'Enable', 'on');
        else
            handles.Current.Tracks1 = handles.Tracking1.Tracks;
        end
        
        
        axes(handles.Channel1data);
        plotTracks(handles.Current.Tracks1, handles.Current.Ind);
        set(handles.showTracksCheck, 'Value', 1);
        set(handles.showTracksCheck, 'Enable', 'on');
        %         set(handles.checkButton,'Enable','on');
        %         set(handles.Track_Preprocess,'Enable','on');
        %         set(handles.HistD_Analysis,'Enable','on');
        %         set(handles.Bound_Analysis,'Enable','on');
        %         set(handles.Export_Params,'Enable','on');
        
    end
    
    
    
    
    
    %store tracking parameters;
    if isfield(Parameters, 'Used');
        if isfield(Parameters.Used,'Tracking')
            handles.Parameters.Used1.Tracking = Parameters.Used.Tracking;
            %         set(handles.lpEdit, 'String', num2str(handles.Parameters.Used.Tracking(1)));
            %         set(handles.hpEdit, 'String', num2str(handles.Parameters.Used.Tracking(2)));
            %         set(handles.thresholdEdit, 'String', num2str(handles.Parameters.Used.Tracking(3)));
            %         set(handles.windowEdit, 'String', num2str(handles.Parameters.Used.Tracking(4)));
            %         if length(handles.Parameters.Used1.Tracking) > 4
            %             set(handles.jumpEdit, 'String', num2str(handles.Parameters.Used.Tracking(5)));
            %             set(handles.gapsCheck, 'Value', handles.Parameters.Used.Tracking(6));
            %             set(handles.shTrackEdit, 'String', num2str(handles.Parameters.Used.Tracking(7)));
            %         end
        end
    end
    
    % check if there are acquisition parameters.
    
    if isfield(Parameters,'Acquisition');
        handles.Parameters.Acquisition1 = Parameters.Acquisition;
        handles.Parameters.Used1.Acquisition = Parameters.Acquisition;
        if isfield(Parameters.Acquisition,'frameTime')
            set(handles.FrameIntervalText,'String',num2str(Parameters.Acquisition.frameTime,3));
        else
            handles.Parameters.Acquisition1.frameTime = 0;
            warndlg('It appears that this dataset was not preprocessed.', 'No Preprocess Data');
            
        end
        if isfield(Parameters.Acquisition,'pixelSize')
            set(handles.PxSizeText,'String',num2str(Parameters.Acquisition.pixelSize,3));
        else
            handles.Parameters.Acquisition1.pixelSize = 0;
        end
    else
        handles.Parameters.Acquisition1 = [];
        warndlg('It appears that this dataset was not preprocessed.', 'No Preprocess Data');
    end
    
    % check if there are checking parameters.
    
    if isfield(Parameters, 'discardCheck');
        handles.Parameters.discardCheck1 = Parameters.discardCheck;
        handles.Parameters.Used1.discardCheck = Parameters.discardCheck;
    end
    
    % check if there are analysis parameters.
    
    if isfield(Parameters, 'Analysis');
        handles.Parameters.Analysis1 = Parameters.Analysis;
        handles.Parameters.Used1.Analysis = Parameters.Analysis;
    end
    
    % check if there are parameters for the bound molecules.
    
    if isfield(Parameters, 'Bound');
        handles.Parameters.Bound1 = Parameters.Bound;
        handles.Parameters.Used1.Bound = Parameters.Bound;
    end
    
    set(handles.LoadData2Push,'Enable','on');
    
    
    
    if isfield(handles,'Data2') && ~isempty(handles.Data2)
        set(handles.ColocPush,'Enable','on');
    end
    handles.isAlign = 0;
end
% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in LoadData2Push.
function LoadData2Push_Callback(hObject, eventdata, handles)
% hObject    handle to LoadData2Push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Choose the default folder
if isfield(handles.Data1,'pathName') && ~isempty(handles.Data1.pathName)
    defFold = handles.Data1.pathName;
else
    defFold = pwd;
end

[FileName,PathName,FilterIndex] = ...
    uigetfile('.mat','Select Matlab File with tracked data',defFold);

if FileName ~=0;
    ListVar = whos('-file',[PathName,FileName]);
else
    return
end
% Reset the data content of the GUI
handles.Process2 = [];
handles.Parameters.Used2 = [];
handles.Tracking2 = [];
handles.Analysis2 = [];
handles.PreAnalysis2 = [];

handles.isFitPSF2 = 0;

IN = load([PathName,FileName],'Version');

if IN.Version == 2.0;
    
    %Load all the contents of the structure.
    IN = load([PathName,FileName], 'Results');
    
    handles.Data2 =  IN.Results.Data;
    
    % Make sure there are the same number of images in Channel 2 as
    % Channel 1
    if handles.Data2.nImages ~= handles.Data1.nImages || ...
            handles.Data2.imageStack(1).height ~= handles.Data1.imageStack(1).height || ...
            handles.Data2.imageStack(1).width ~= handles.Data1.imageStack(1).width
        msgbox('Size of Channel 2 images differs from Channel 1. Please choose a different file');
        return
    end;
    
    
    
%     handles.Data2.fileName = IN.Results.Data.fileName;
%     handles.Data2.pathName = IN.Results.Data.pathName;
    Parameters = IN.Results.Parameters;
    handles.Process2 = IN.Results.Process;
    handles.Tracking2 = IN.Results.Tracking;
    handles.isFitPSF2 = IN.Results.isFitPSF;
    handles.PreAnalysis2 = IN.Results.PreAnalysis;
    
    % Disable some of the controls
    % [They will be re-enabled if the variables are found within the input file]
    
    %     % Disable the otpions to visualize hand checked tracks
    %     set(handles.TrackViewPop, 'Value',1);
    %     set(handles.TrackViewPop, 'Enable', 'off');
    %
    %     % Disable the options to visualize filtered images
    %     set(handles.ImageViewPop, 'Value',1);
    %     set(handles.ImageViewPop, 'Enable', 'off');
    %
    %     % Disable the options to visualize localized particles
    %     set(handles.showParticlesCheck, 'Value', 0);
    %     set(handles.showParticlesCheck, 'Enable', 'off');
    %
    %     % Disable the options to visualize the tracks
    %     set(handles.showTracksCheck, 'Value', 0);
    %     set(handles.showTracksCheck, 'Enable', 'off');
    
    %Enable filter button & roiButton
    %     set(handles.filterButton,'Enable','on');
    %     set(handles.roiButton,'Enable','on');
    %     set(handles.lpEdit,'Enable','on');
    %     set(handles.hpEdit,'Enable','on');
    
    
    %     % reset current index.
    %     handles.Current = [];
    %     handles.Current.Ind = 1;
    
    % Select current image and set colormap limits.
    %calculate the minimum & maximum values in the stack
    Mn_o = 100000;
    Mx_o = 0;
    
    for i = 1:handles.Data2.nImages
        if min(handles.Data2.imageStack(i).data(:)) < Mn_o
            Mn_o = min(handles.Data2.imageStack(i).data(:));
        end
        if max(handles.Data2.imageStack(i).data(:)) > Mx_o
            Mx_o = max(handles.Data2.imageStack(i).data(:));
        end
    end
    %
    handles.Current.clims2 = [Mn_o Mx_o];
    handles.Data2.clims_orig = handles.Current.clims2;
    set(handles.BlackSlider_C2,'Min',0,'Max',2*double(Mx_o),'SliderStep',[1/((2*double(Mx_o))-1),10/((2*double(Mx_o))-1)],'Value',double(Mn_o));
    set(handles.WhiteSlider_C2,'Min',0,'Max',2*double(Mx_o),'SliderStep',[1/((2*double(Mx_o))-1),10/((2*double(Mx_o))-1)],'Value',double(Mx_o));
    set(handles.BlackEdit_C2,'String',num2str(Mn_o));
    set(handles.WhiteEdit_C2,'String',num2str(Mx_o));
    
    %Enable other fields & plot if the data is available
    if isfield(handles.Process2,'filterStack')
        
        Mn_f = 100000;
        Mx_f = 0;
        for i = 1:handles.Data2.nImages
            if min(handles.Process2.filterStack(i).data(:)) < Mn_f
                Mn_f = min(handles.Process2.filterStack(i).data(:));
            end
            if max(handles.Process2.filterStack(i).data(:)) > Mx_f
                Mx_f = max(handles.Process2.filterStack(i).data(:));
            end
            %             if Mx_f/Mn_f > 6
            %                 handles.Data2.clims_filt = [Mn_f 0.25*Mx_f];
            %             elseif Mx_o/Mn_o > 4
            %                 handles.Data2.clims_filt = [Mn_f 0.5*Mx_f];
            %             else
            handles.Data2.clims_filt = [Mn_f Mx_f];
            %             end
            
        end
        
    end
    %%%Add check on ImageViewPop
    if get(handles.ImageViewPop,'Value') == 2 && isfield(handles.Process2,'filterStack')
        handles.Current.Image2 = handles.Process2.filterStack(handles.Current.Ind).data;
        handles.Current.clims2 = handles.Data2.clims_filt;
    else
        handles.Current.Image2 = handles.Data2.imageStack(handles.Current.Ind).data;
        handles.Current.clims2 = handles.Data2.clims_orig;
    end
    
    
    % Display current Image
    
    
    axes(handles.Channel2data);
    colormap(gray);
    imagesc(handles.Current.Image2, handles.Current.clims2);
    axis image;
    
    
    
    
    
    % Check if there is the filtered stack
    if isfield(handles.Process2, 'filterStack') && strcmpi(get(handles.ImageViewPop,'Enable'),'off')
        % Enable the options to visualize filtered images
        set(handles.ImageViewPop, 'Value',1);
        set(handles.ImageViewPop, 'Enable', 'on');
        %         set(handles.thresholdEdit,'Enable','on');
        %         set(handles.windowEdit,'Enable','on');
        %         set(handles.fitPSF,'Enable','on');
        %         set(handles.findpeakButton,'Enable','on');
        
        
    end
    %Check if there is ROI data
    if isfield(handles.Process2,'ROIpos')
        nROIs2 = length(handles.Process2.ROIpos);
        
        if get(handles.showROIsCheck,'Value');
            axes(handles.Channel2data);
            plotROI(handles.Process2.ROIpos);
        end
        
        %         set(handles.roiRemovePush,'Enable','on');
        %         set(handles.roiNamePush,'Enable','on');
    end
    
    
    
    % Check if there is the data about the particles;
    if isfield(handles.Tracking2, 'Centroids')||isfield(handles.Tracking1,'Particles')
        % Check if there particles are fit with gaussians
        
        
        if isfield(handles.Tracking2, 'Particles');
            
            Centroids =  handles.Tracking2.Particles;
            
            
        else
            
            Centroids =  handles.Tracking1.Centroids;
            
        end
        
        
        if get(handles.showParticlesCheck,'Value');
            axes(handles.Channel2data);
            plotParticle(Centroids,handles.Current.Ind,handles.isFitPSF2);
        end
        %         handles.arePeaksFiltered = 0;
        %         set(handles.Intensity_thr,'Enable','on');
        %         set(handles.sigmaLow,'Enable','on');
        %         set(handles.sigmaHigh,'Enable','on');
        %         set(handles.filterPeaks,'Enable','on');
        %         set(handles.editPeakButton,'Enable','on');
        %         set(handles.xPlot,'Enable','on');
        %         set(handles.yPlot,'Enable','on');
        %         set(handles.PlotVariables,'Enable','on');
        %         set(handles.jumpEdit,'Enable','on');
        %         set(handles.sigmaLow,'Enable','on');
        %         set(handles.shTrackEdit,'Enable','on');
        %         set(handles.gapsCheck,'Enable','on');
        %         set(handles.trackButton,'Enable','on');
        
    end
    
    % Check if there is information about tracks;
    if isfield(handles.Tracking2, 'Tracks')
        % Check if there are hand checked tracks
        if isfield(handles.Tracking2, 'CheckTracks') && get(handles.TrackViewPop,'Value') == 2
            handles.Current.Tracks2 = handles.Tracking2.CheckTracks;
            
        else
            handles.Current.Tracks2 = handles.Tracking2.Tracks;
        end
        
        
        
        if get(handles.showTracksCheck, 'Value');
            
            axes(handles.Channel2data);
            plotTracks(handles.Current.Tracks2, handles.Current.Ind);
        end
        %         set(handles.checkButton,'Enable','on');
        %         set(handles.Track_Preprocess,'Enable','on');
        %         set(handles.HistD_Analysis,'Enable','on');
        %         set(handles.Bound_Analysis,'Enable','on');
        %         set(handles.Export_Params,'Enable','on');
        
    end
    
    
    
    
    
    %store tracking parameters;
    if isfield(Parameters, 'Tracking');
        handles.Parameters.Used2.Tracking = Parameters.Tracking;
        %         set(handles.lpEdit, 'String', num2str(handles.Parameters.Used.Tracking(1)));
        %         set(handles.hpEdit, 'String', num2str(handles.Parameters.Used.Tracking(2)));
        %         set(handles.thresholdEdit, 'String', num2str(handles.Parameters.Used.Tracking(3)));
        %         set(handles.windowEdit, 'String', num2str(handles.Parameters.Used.Tracking(4)));
        %         if length(handles.Parameters.Used1.Tracking) > 4
        %             set(handles.jumpEdit, 'String', num2str(handles.Parameters.Used.Tracking(5)));
        %             set(handles.gapsCheck, 'Value', handles.Parameters.Used.Tracking(6));
        %             set(handles.shTrackEdit, 'String', num2str(handles.Parameters.Used.Tracking(7)));
        %         end
    end
    
    % check if there are acquisition parameters.
    
    if isfield(Parameters,'Acquisition');
        handles.Parameters.Acquisition2 = Parameters.Acquisition;
        handles.Parameters.Used2.Acquisition = Parameters.Acquisition;
        if isfield(Parameters.Acquisition,'frameTime')
            if isfield(handles.Parameters.Acquisition1,'frameTime') && handles.Parameters.Acquisition2.frameTime ~= handles.Parameters.Acquisition1.frameTime
                ChanInterQuest = questdlg('Frame Intervals are different for the 2 channels, which one do you want to use?','Channel 1','Channel 2','Channel 1');
                if strcmp(ChanInterQuest, 'Channel 2')
                    set(handles.FrameIntervalText,'String',num2str(Parameters.Acquisition.frameTime,3));
                end
            else
                set(handles.FrameIntervalText,'String',num2str(Parameters.Acquisition.frameTime,3));
            end
        else
            warndlg('It appears that this dataset was not preprocessed.', 'No Preprocess Data');
        end
        
    end
    
    % check if there are checking parameters.
    
    if isfield(Parameters, 'discardCheck');
        handles.Parameters.discardCheck2 = Parameters.discardCheck;
        handles.Parameters.Used2.discardCheck = Parameters.discardCheck;
    end
    
    % check if there are analysis parameters.
    
    if isfield(Parameters, 'Analysis');
        handles.Parameters.Analysis2 = Parameters.Analysis;
        handles.Parameters.Used2.Analysis = Parameters.Analysis;
    end
    
    % check if there are parameters for the bound molecules.
    
    if isfield(Parameters, 'Bound');
        handles.Parameters.Bound2 = Parameters.Bound;
        handles.Parameters.Used2.Bound = Parameters.Bound;
    end
    
    %Create a merged image from the 2 channels with Channel 1 red & Channel
    %2 green,
    handles.Data2.imageStack(handles.Current.Ind).data;
    handles.Current.mergedImage = merge2colorData(handles.Current.Image1,handles.Current.Image2,[1 2],handles.Current.clims1,handles.Current.clims2);
    handles.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,1),size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages);
    handles.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,1),size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages);
    for i = 1:handles.Data1.nImages
        handles.mergedStack.Image1(:,:,i) = handles.Data1.imageStack(i).data;
        handles.mergedStack.Image2(:,:,i) = handles.Data2.imageStack(i).data;
    end
    handles.Current.mergedStack = handles.mergedStack;
    handles.Current.MergedInd = handles.Current.Ind;
    set(handles.MergedData,'Visible','on');
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    set(handles.MergeChannel1Check,'Visible','on');
    set(handles.MergeChannel2Check,'Visible','on');
    set(handles.MergeChannel1Check,'Value',1);
    set(handles.MergeChannel2Check,'Value',1);
    set(handles.ViewSelectStatic,'Visible','on');
    set(handles.MergeViewPop,'Visible','on');
    set(handles.MergeViewPop,'Enable','on');
    set(handles.MergeViewPop,'Value',1);
    set(handles.ExportMergeImage,'Enable','on');
    
    axis image;
    set(handles.mergeSlider,'Visible','off');
    
    set(handles.ChanAlignPush,'Enable','on');
    if isfield(handles,'Data1') && ~isempty(handles.Data1)
        set(handles.ColocPush,'Enable','on');
    end
    handles.isAlign = 0;
    ax(1) = handles.Channel1data;
    ax(2) = handles.Channel2data';
    ax(3) = handles.MergedData';
    linkaxes(ax);
    
    % Update handles structure
    guidata(hObject,handles);
    
end

% --- Executes on selection change in ImageViewPop.
function ImageViewPop_Callback(hObject, eventdata, handles)
% hObject    handle to ImageViewPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ImageViewPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ImageViewPop
ax = [handles.Channel1data, handles.Channel2data, handles.MergedData];
linkaxes(ax,'off');

whatToDraw = get(hObject, 'Value');
% Image the selected image
switch whatToDraw
    case 1
        handles.Current.Image1 = ...
            handles.Data1.imageStack(handles.Current.Ind).data;
        
        handles.Current.clims1 = handles.Data1.clims_orig;
        if isfield(handles,'Data2')
            handles.Current.Image2 = ...
                handles.Data2.imageStack(handles.Current.Ind).data;
            handles.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,1),size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages);
            handles.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,1),size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages);
            for i = 1:handles.Data1.nImages
                handles.mergedStack.Image1(:,:,i) = handles.Data1.imageStack(i).data;
                handles.mergedStack.Image2(:,:,i) = handles.Data2.imageStack(i).data;
            end
            handles.Current.clims2 = handles.Data2.clims_orig;
        end
        
    case 2
        handles.Current.Image1 = ...
            handles.Process1.filterStack(handles.Current.Ind).data;
        handles.Current.clims1 = handles.Data1.clims_filt;
        
        if isfield(handles,'Process2')
            handles.Current.Image2 = ...
                handles.Process2.filterStack(handles.Current.Ind).data;
            handles.mergedStack.Image1 = zeros(size(handles.Process1.filterStack(1).data,1),size(handles.Process1.filterStack(1).data,2),handles.Data1.nImages);
            handles.mergedStack.Image2 = zeros(size(handles.Process2.filterStack(1).data,1),size(handles.Process2.filterStack(1).data,2),handles.Data2.nImages);
            for i = 1:handles.Data1.nImages
                handles.mergedStack.Image1(:,:,i) = handles.Process1.filterStack(i).data;
                handles.mergedStack.Image2(:,:,i) = handles.Process2.filterStack(i).data;
            end
            handles.Current.clims2 = handles.Data2.clims_filt;
        end
end
set(handles.BlackSlider_C1,'Min',0.0,'Max',2*double(handles.Current.clims1(2)),'SliderStep',[1/(2*double(handles.Current.clims1(2)) - 1), 10/(2*double(handles.Current.clims1(2))-1)]);
set(handles.BlackSlider_C1,'Val',handles.Current.clims1(1));
set(handles.BlackEdit_C1,'String',num2str(handles.Current.clims1(1)));

set(handles.WhiteSlider_C1,'Min',0.0,'Max',2*double(handles.Current.clims1(2)),'SliderStep',[1/(2*double(handles.Current.clims1(2)) - 1), 10/(2*double(handles.Current.clims1(2))-1)]);
set(handles.WhiteSlider_C1,'Val',handles.Current.clims1(2));
set(handles.WhiteEdit_C1,'String',num2str(handles.Current.clims1(2)));

if isfield(handles,'Process2')
    set(handles.BlackSlider_C2,'Min',0.0,'Max',2*double(handles.Current.clims2(2)),'SliderStep',[1/(2*double(handles.Current.clims2(2)) - 1), 10/(2*double(handles.Current.clims2(2))-1)]);
    set(handles.BlackSlider_C2,'Val',handles.Current.clims2(1));
    set(handles.BlackEdit_C2,'String',num2str(handles.Current.clims2(1)));
    
    set(handles.WhiteSlider_C2,'Min',0.0,'Max',2*double(handles.Current.clims2(2)),'SliderStep',[1/(2*double(handles.Current.clims2(2)) - 1), 10/(2*double(handles.Current.clims2(2))-1)]);
    set(handles.WhiteSlider_C2,'Val',handles.Current.clims2(2));
    set(handles.WhiteEdit_C2,'String',num2str(handles.Current.clims2(2)));
end

% handles.Current.clims1 = [min(handles.Current.Image1(:)) max(handles.Current.Image1(:))];
axes(handles.Channel1data);
imagesc(handles.Current.Image1, handles.Current.clims1);
axis image;

if isfield(handles.Current,'Image2')
    %     handles.Current.clims2 = [min(handles.Current.Image2(:)) max(handles.Current.Image2(:))];
    axes(handles.Channel2data);
    imagesc(handles.Current.Image2, handles.Current.clims2);
    axis image;
    
    switch get(handles.MergeViewPop,'Value')
        case 1
            handles.Current.mergedStack.Image1 = handles.mergedStack.Image1;
            handles.Current.mergedStack.Image2 = handles.mergedStack.Image2;
            handles.Current.MergedInd = handles.Current.Ind;
        case 2
            handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[2 3 1]);
            handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[2 3 1]);
        case 3
            handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[1 3 2]);
            handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[1 3 2]);
    end
    handles.Current.mergedImage = merge2colorData(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd),handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd),[1 2],handles.Current.clims1,handles.Current.clims2);
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    if get(handles.MergeViewPop,'Value') == 1
        axis image
    else
        ax(3) = [];
        
    end
end
linkaxes(ax);
if get(handles.showParticlesCheck, 'Value')
    axes(handles.Channel1data);
    if handles.isFitPSF1
        
        plotParticle(handles.Tracking1.Particles, ...
            handles.Current.Ind, handles.isFitPSF1);
    else
        plotParticle(handles.Tracking1.Centroids, ...
            handles.Current.Ind, handles.isFitPSF1);
    end
    
    if isfield(handles,'Tracking2')
        axes(handles.Channel2data);
        if handles.isFitPSF2
            
            plotParticle(handles.Tracking2.Particles, ...
                handles.Current.Ind, handles.isFitPSF2);
        else
            plotParticle(handles.Tracking2.Centroids, ...
                handles.Current.Ind, handles.isFitPSF2);
        end
    end
end

if get(handles.showTracksCheck, 'Value')
    axes(handles.Channel1data);
    plotTracks(handles.Current.Tracks1,handles.Current.Ind);
    
    
    if isfield(handles.Current,'Tracks2')
        axes(handles.Channel2data);
        plotTracks(handles.Current.Tracks2,handles.Current.Ind);
    end
end

if isfield(handles.Process1,'ROIpos') && get(handles.showROIsCheck,'Value')
    axes(handles.Channel1data);
    plotROI(handles.Process1.ROIpos);
    
end

if isfield(handles,'Process2')
    if isfield(handles.Process2,'ROIpos') && get(handles.showROIsCheck,'Value')
        axes(handles.Channel2data);
        plotROI(handles.Process2.ROIpos);
    end
end

% Update handles structure
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function ImageViewPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageViewPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in TrackViewPop.
function TrackViewPop_Callback(hObject, eventdata, handles)
% hObject    handle to TrackViewPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns TrackViewPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TrackViewPop
whatToDraw = get(hObject, 'Value');
% Image the selected image
switch whatToDraw
    case 1
        handles.Current.Tracks1 = handles.Tracking1.Tracks;
        if isfield(handles,'Tracking2')
            handles.Current.Tracks2 = handles.Tracking2.Tracks;
        end
        
    case 2
        handles.Current.Tracks1 = handles.Tracking1.CheckTracks;
        if isfield(handles,'Tracking2')
            handles.Current.Tracks2 = handles.Tracking2.CheckTracks;
            
        end
end;
axes(handles.Channel1data);
plotTracks(handles.Current.Tracks1,handles.Current.Ind);
if isfield(handles.Current,'Tracks2')
    axes(handles.Channel2data);
    plotTracks(handles.Current.Tracks2,handles.Current.Ind);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function TrackViewPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrackViewPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function imageSlider_Callback(hObject, eventdata, handles)
% hObject    handle to imageSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
ax = [handles.Channel1data, handles.Channel2data, handles.MergedData];
linkaxes(ax,'off');

whatToDraw = get(handles.ImageViewPop, 'Value');
% Image the selected image
handles.Current.Ind = round(get(hObject,'Value'));
switch whatToDraw
    case 1
        handles.Current.Image1 = ...
            handles.Data1.imageStack(handles.Current.Ind).data;
        axes(handles.Channel1data);
        imagesc(handles.Current.Image1, handles.Current.clims1);
        axis image;
        if isfield(handles,'Data2')
            handles.Current.Image2 = handles.Data2.imageStack(handles.Current.Ind).data;
            axes(handles.Channel2data);
            imagesc(handles.Current.Image2, handles.Current.clims2);
            axis image;
            
            %display the merged image
            if get(handles.MergeViewPop,'Value') == 1
                if get(handles.MergeChannel1Check,'Value')
                    mergeChannel1 = handles.Current.Image1;
                else
                    mergeChannel1 = zeros(size(handles.Current.Image1));
                end
                if get(handles.MergeChannel2Check,'Value')
                    mergeChannel2 = handles.Current.Image2;
                else
                    mergeChannel2 = zeros(size(handles.Current.Image2));
                end
                handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
                axes(handles.MergedData);
                image(handles.Current.mergedImage);
                axis image
            end
        end
    case 2
        handles.Current.Image1 = ...
            handles.Process1.filterStack(handles.Current.Ind).data;
        axes(handles.Channel1data);
        ax(1) = gca;
        imagesc(handles.Current.Image1,handles.Current.clims1);
        axis image;
        if isfield(handles,'Process2')
            handles.Current.Image2 = handles.Process2.filterStack(handles.Current.Ind).data;
            axes(handles.Channel2data);
            ax(2) = gca;
            imagesc(handles.Current.Image2, handles.Current.clims2);
            axis image;
            if get(handles.MergeViewPop,'Value') == 1
                if get(handles.MergeChannel1Check,'Value')
                    mergeChannel1 = handles.Current.Image1;
                else
                    mergeChannel1 = zeros(size(handles.Current.Image1));
                end
                if get(handles.MergeChannel2Check,'Value')
                    mergeChannel2 = handles.Current.Image2;
                else
                    mergeChannel2 = zeros(size(handles.Current.Image2));
                end
                handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
                axes(handles.MergedData);
                ax(3) = gca;
                image(handles.Current.mergedImage);
                axis image
            end
            
        end
end;



%Plot particles for channel 1 if requested & channel 2 if available
if get(handles.showParticlesCheck, 'Value')
    if handles.isFitPSF1
        axes(handles.Channel1data);
        plotParticle(handles.Tracking1.Particles, ...
            handles.Current.Ind, handles.isFitPSF1);
        
    else
        axes(handles.Channel1data);
        plotParticle(handles.Tracking1.Centroids, ...
            handles.Current.Ind, handles.isFitPSF1);
    end
    
    if isfield(handles,'Tracking2') && isfield(handles,'isFitPSF2')
        if handles.isFitPSF2
            axes(handles.Channel2data);
            plotParticle(handles.Tracking2.Particles, ...
                handles.Current.Ind, handles.isFitPSF2);
            
        else
            axes(handles.Channel2data);
            plotParticle(handles.Tracking2.Centroids, ...
                handles.Current.Ind, handles.isFitPSF2);
        end
    end
end

%Plot tracks if requested & available
if get(handles.showTracksCheck, 'Value')
    axes(handles.Channel1data);
    plotTracks(handles.Current.Tracks1, handles.Current.Ind);
    if isfield(handles.Current,'Tracks2')
        axes(handles.Channel2data);
        plotTracks(handles.Current.Tracks2, handles.Current.Ind);
    end
end

if get(handles.showROIsCheck,'Value')
    if isfield(handles.Process1,'ROIpos')
        axes(handles.Channel1data);
        plotROI(handles.Process1.ROIpos);
    end
    if isfield(handles,'Process2');
        if isfield(handles.Process2,'ROIpos');
            axes(handles.Channel2data);
            plotROI(handles.Process2.ROIpos);
        end
    end
    
    
end

% Update the image indicator in the panel
set(handles.ImageCount, 'String', ['Image ', num2str(handles.Current.Ind),'/',num2str(handles.Data1.nImages)]);

if get(handles.LocationDrawCheck,'Value')
    if get(handles.MergeViewPop,'Value') == 2
        tyline = [0 size(handles.Current.Image1,2)];
        txline = [handles.Current.Ind handles.Current.Ind];
        xline = [0 size(handles.Current.Image1,2)];
        yline = [handles.Current.MergedInd handles.Current.MergedInd];
    elseif get(handles.MergeViewPop,'Value') == 3
        tyline = [0 size(handles.Current.Image1,1)];
        txline = [handles.Current.Ind handles.Current.Ind];
        yline = [0 size(handles.Current.Image1,1)];
        xline = [handles.Current.MergedInd handles.Current.MergedInd];
    end
    axes(handles.Channel1data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    axes(handles.Channel2data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    
    axes(handles.MergedData);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(txline,tyline,'w','LineWidth',0.75);
    hold off
end
linkaxes(ax);

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


% --- Executes on selection change in MergeViewPop.
function MergeViewPop_Callback(hObject, eventdata, handles)
% hObject    handle to MergeViewPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MergeViewPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MergeViewPop
ax = [handles.Channel1data, handles.Channel2data, handles.MergedData];
linkaxes(ax,'off');

showWhich = get(hObject,'Value');
turnOnOff1 = get(handles.MergeChannel1Check,'Value');
turnOnOff2 = get(handles.MergeChannel2Check,'Value');

switch showWhich
    case 1 %X*Y
        set(handles.mergeSlider,'Visible','off');
        %         handles.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,1),size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages);
        %         handles.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,1),size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages);
        %         for i = 1:handles.Data1.nImages
        %             if get(handles.ImageViewPop,'Value') == 1
        %                 handles.mergedStack.Image1(:,:,i) = handles.Data1.imageStack(i).data;
        %                 handles.mergedStack.Image2(:,:,i) = handles.Data2.imageStack(i).data;
        %             else
        %                 handles.mergedStack.Image1(:,:,i) = handles.Process1.filterStack(i).data;
        %                 handles.mergedStack.Image2(:,:,i) = handles.Process2.filterStack(i).data;
        %             end
        %         end
        handles.Current.mergedStack = handles.mergedStack;
        handles.Current.MergedInd = handles.Current.Ind;
        if turnOnOff1
            mergeChannel1 = handles.Current.Image1;
        else
            mergeChannel1 = zeros(size(handles.Current.Image1));
        end
        if turnOnOff2
            mergeChannel2 = handles.Current.Image2;
        else
            mergeChannel2 = zeros(size(handles.Current.Image2));
        end
        handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
        axes(handles.MergedData);
        
        image(handles.Current.mergedImage);
        axis image
        set(handles.SliceCount,'Visible','off');
        set(handles.LocationDrawCheck,'Visible','off');
        set(handles.LocationDrawCheck,'Value',0);
        linkaxes(ax);
        
    case 2 %X*T
        set(handles.mergeSlider,'Visible','on');
        set(handles.mergeSlider,'Enable','on');
        
        set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,1));
        
        set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,1));
        set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,1)-1) 1/(size(handles.Data1.imageStack(1).data,1)-1)]);
        handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[2 3 1]);
        handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[2 3 1]);
        
        %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,1));
        %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,1));
        %         for i = 1:size(handles.Data2.imageStack(1).data,1)
        %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(i,:,:),size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages);
        %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(i,:,:),size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages);
        %         end
        handles.Current.MergedInd = 1;
        if turnOnOff1
            mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
        else
            mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
        end
        if turnOnOff2
            mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
        else
            mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
        end
        handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
        axes(handles.MergedData);
        image(handles.Current.mergedImage);
        xlabel('Frame Number');
        %         axis image
        set(handles.SliceCount,'Visible','on');
        set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,1))]);
        set(handles.LocationDrawCheck,'Visible','on');
        ax(3) = [];
        linkaxes(ax);
        
    case 3 %Y*T
        set(handles.mergeSlider,'Visible','on');
        set(handles.mergeSlider,'Enable','on');
        
        set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,2))
        set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,2));
        set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,2)-1) 1/(size(handles.Data1.imageStack(1).data,2)-1)]);
        handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[1 3 2]);
        handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[1 3 2]);
        %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,2));
        %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,2));
        %         for i = 1:size(handles.Data2.imageStack(1).data,2)
        %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(:,i,:),size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages);
        %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(:,i,:),size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages);
        %         end
        handles.Current.MergedInd = 1;
        if turnOnOff1
            mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
        else
            mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
        end
        if turnOnOff2
            mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
        else
            mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
        end
        handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
        axes(handles.MergedData);
        image(handles.Current.mergedImage);
        xlabel('Frame Number');
        %         axis image
        set(handles.SliceCount,'Visible','on');
        set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,2))]);
        set(handles.LocationDrawCheck,'Visible','on');
        ax(3) = [];
        linkaxes(ax);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function MergeViewPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MergeViewPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in MergeChannel1Check.
function MergeChannel1Check_Callback(hObject, eventdata, handles)
% hObject    handle to MergeChannel1Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MergeChannel1Check
turnOnOff1 = get(hObject,'Value');
turnOnOff2 = get(handles.MergeChannel2Check,'Value');

if turnOnOff1
    mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
else
    mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
end
if turnOnOff2
    mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
else
    mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
end
handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
axes(handles.MergedData);
image(handles.Current.mergedImage);

if get(handles.MergeViewPop, 'Value') > 1
    xlabel('Frame Number');
else
    axis image
end
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in MergeChannel2Check.
function MergeChannel2Check_Callback(hObject, eventdata, handles)
% hObject    handle to MergeChannel2Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MergeChannel2Check

ax = [handles.Channel1data, handles.Channel2data, handles.MergedData];
linkaxes(ax,'off');

turnOnOff2 = get(hObject,'Value');
turnOnOff1 = get(handles.MergeChannel1Check,'Value');

if turnOnOff1
    mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
else
    mergeChannel1 = zeros(size(handles.mergedStack.Current.Image1(:,:,handles.Current.MergedInd)));
end
if turnOnOff2
    mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
else
    mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
end
handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
axes(handles.MergedData);
image(handles.Current.mergedImage);

if get(handles.MergeViewPop, 'Value') > 1
    ax(3) = [];
    xlabel('Frame Number');
else
    axis image
end
linkaxes(ax);

% Update handles structure
guidata(hObject, handles);


% --- Executes on slider movement.
function mergeSlider_Callback(hObject, eventdata, handles)
% hObject    handle to mergeSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if get(handles.MergeViewPop,'Value') == 2
    handles.Current.MergedInd = (size(handles.Data1.imageStack(1).data,1)+1) - round(get(hObject,'Value'));
else
    handles.Current.MergedInd = (size(handles.Data1.imageStack(1).data,2)+1) - round(get(hObject,'Value'));
end
turnOnOff1 = get(handles.MergeChannel1Check,'Value');
turnOnOff2 = get(handles.MergeChannel2Check,'Value');
if turnOnOff1
    mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
else
    mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
end
if turnOnOff2
    mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
else
    mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
end
handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
axes(handles.MergedData);
image(handles.Current.mergedImage);
if get(handles.MergeViewPop, 'Value') > 1
    xlabel('Frame Number');
end
% axis image
set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.Current.mergedStack.Image1,3))]);
if get(handles.LocationDrawCheck,'Value')
    if get(handles.MergeViewPop,'Value') == 2
        xline = [0 size(handles.Current.Image1,2)];
        yline = [handles.Current.MergedInd handles.Current.MergedInd];
        tyline = xline;
    elseif get(handles.MergeViewPop,'Value') == 3
        yline = [0 size(handles.Current.Image1,1)];
        xline = [handles.Current.MergedInd handles.Current.MergedInd];
        tyline = yline;
    end
    txline = [handles.Current.Ind handles.Current.Ind];
    axes(handles.Channel1data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    
    axes(handles.Channel2data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    
    axes(handles.MergedData)
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(txline,tyline,'w','LineWidth',0.75);
    hold off
end
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function mergeSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mergeSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in LocationDrawCheck.
function LocationDrawCheck_Callback(hObject, eventdata, handles)
% hObject    handle to LocationDrawCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LocationDrawCheck
if get(hObject,'Value')
    if get(handles.MergeViewPop,'Value') == 2
        xline = [0 size(handles.Current.Image1,2)];
        yline = [handles.Current.MergedInd handles.Current.MergedInd];
        tyline = xline;
    elseif get(handles.MergeViewPop,'Value') == 3
        yline = [0 size(handles.Current.Image1,1)];
        xline = [handles.Current.MergedInd handles.Current.MergedInd];
        tyline = yline;
    end
    txline = [handles.Current.Ind handles.Current.Ind];
    axes(handles.Channel1data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    
    axes(handles.Channel2data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    
    axes(handles.MergedData)
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(txline,tyline,'w','LineWidth',0.75);
    hold off
else
    axes(handles.Channel1data);
    delete(findobj(gca,'LineWidth',0.75));
    
    
    axes(handles.Channel2data);
    delete(findobj(gca,'LineWidth',0.75));
    
    
    axes(handles.MergedData)
    delete(findobj(gca,'LineWidth',0.75));
    
end


% --------------------------------------------------------------------
function File_menu_Callback(hObject, eventdata, handles)
% hObject    handle to File_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Edit_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ChanColor_Callback(hObject, eventdata, handles)
% hObject    handle to ChanColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Save_Mat_Callback(hObject, eventdata, handles)
% hObject    handle to Save_Mat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Assign .mat file where variables will be saved.
defName = [handles.Data1.pathName, handles.Data1.fileName];
SearchStr = '(.*)\.\w*';
defName = regexprep(defName, SearchStr, '$1');
FilterSpec = {'*.mat'};
[FileNameOut,PathNameOut] = uiputfile(FilterSpec,'Save Variables',defName);
if FileNameOut ~=0
    Version = 1.0;
    
    %     save([PathNameOut, FileNameOut], 'Version');
    Results.Data1 = handles.Data1;
    Results.Data2 = handles.Data2;
    %     Results.Current = handles.Current;
    Results.Parameters = handles.Parameters;
    Results.isAlign = handles.isAlign;
    Results.Process1 = handles.Process1;
    Results.Process2 = handles.Process2;
    Results.Tracking1 = handles.Tracking1;
    Results.Tracking2 = handles.Tracking2;
    %     Results.Analysis1 = handles.Analysis1;
    %     Results.Analysis2 = handles.Analysis2;
    Results.PreAnalysis1 = handles.PreAnalysis1;
    Results.PreAnalysis2 = handles.PreAnalysis2;
    Results.isFitPSF1 = handles.isFitPSF1;
    Results.isFitPSF2 = handles.isFitPSF2;
    %     Results.mergedStack = handles.mergedStack;
    save([PathNameOut, FileNameOut], 'Version','Results', '-v7.3');
    
end

% --------------------------------------------------------------------
function Load_Mat_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Mat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Open MatLab file and get variables
[FileName,PathName,FilterIndex] = ...
    uigetfile('.mat','Select Matlab File with tracked data');

if FileName ~=0;
    ListVar = whos('-file',[PathName,FileName]);
else
    return
end
% Reset the data content of the GUI
handles.Data1 = [];
handles.Data2 = [];
%     Results.Current = handles.Current;
handles.Parameters = [];
handles.isAlign = 0;
handles.Process1 = [];
handles.Process2 = [];
handles.Tracking1 = [];
handles.Tracking2 = [];
handles.PreAnalysis1 = [];
handles.PreAnalysis2 = [];
handles.isFitPSF1 = 0;
handles.isFitPSF2 = 0;

handles.Current = [];
handles.Current.Ind = 1;



if ismember('Version',{ListVar.name});
    % THIS IS IF THE MAT FILE HAS BEEN CREATED WITH THE CURRENT VERSION OF
    % THE SOFTWARE.
    IN = load([PathName,FileName],'Version');
    if IN.Version == 1.0;
        IN = load([PathName,FileName], 'Results');
        Results = IN.Results;
        handles.Data1 = Results.Data1;
        handles.Data2 = Results.Data2;
        handles.Process1 = Results.Process1;
        handles.Process2 = Results.Process2;
        if ~isempty(handles.Data1)
            set(handles.LoadData2Push,'Enable','on');
            if ~isempty(handles.Data2)
                set(handles.ChanAlignPush,'Enable','on');
                set(handles.ColocPush,'Enable','on');
            end
        end
        if ~isfield(handles.Data1,'imageStack')
            handles.Data1.imageStack = TIFread(fullfile(handles.Data1.pathName,handles.Data1.fileName));
            handles.Data2.imageStack = TIFread(fullfile(handles.Data2.pathName,handles.Data2.fileName));
            for i = 1:handles.Data1.nImages
                handles.Process1.filterStack(i).data = ...
                    bpass(handles.Data1.imageStack(i).data, 1, 5);
            end
            for i = 1:handles.Data2.nImages
                handles.Process2.filterStack(i).data = ...
                    bpass(handles.Data2.imageStack(i).data, 1, 5);
            end
        end
        %     Results.Current = handles.Current;
        handles.Parameters = Results.Parameters;
        set(handles.PxSizeText,'String',num2str(handles.Parameters.Acquisition1.pixelSize));
        set(handles.FrameIntervalText,'String',num2str(handles.Parameters.Acquisition1.frameTime));
        handles.isAlign = Results.isAlign;
        
        handles.Tracking1 = Results.Tracking1;
        % Check if there are hand checked tracks
        if isfield(handles.Tracking1, 'CheckTracks')
            handles.Current.Tracks1 = handles.Tracking1.CheckTracks;
            set(handles.TrackViewPop, 'Value',2);
            set(handles.TrackViewPop, 'Enable', 'on');
        else
            handles.Current.Tracks1 = handles.Tracking1.Tracks;
            set(handles.showTracksCheck,'Enable','on');
            set(handles.showTracksCheck,'Value',1);
        end
        if isfield(handles.Tracking1,'Particles')
            set(handles.showParticlesCheck,'Enable','on');
            set(handles.showParticlesCheck,'Value',1);
        end
        if isfield(handles.Process1,'ROIpos')
            set(handles.showROIsCheck,'Enable','on');
            set(handles.showROIsCheck,'Value',1);
        end
        handles.Tracking2 = Results.Tracking2;
        % Check if there are hand checked tracks
        if isfield(handles.Tracking2, 'CheckTracks')
            handles.Current.Tracks2 = handles.Tracking2.CheckTracks;
            set(handles.TrackViewPop, 'Value',2);
            set(handles.TrackViewPop, 'Enable', 'on');
        else
            handles.Current.Tracks2 = handles.Tracking2.Tracks;
        end
        handles.PreAnalysis1 = Results.PreAnalysis1;
        handles.PreAnalysis2 = Results.PreAnalysis2;
        handles.isFitPSF1 = Results.isFitPSF1;
        handles.isFitPSF2 = Results.isFitPSF2;
        for i = 1:handles.Data1.nImages
            handles.mergedStack.Image1(:,:,i) = handles.Data1.imageStack(i).data;
            handles.mergedStack.Image2(:,:,i) = handles.Data2.imageStack(i).data;
        end
        %calculate the minimum & maximum values in the stack
        Mn_o = 100000;
        Mx_o = 0;
        
        for i = 1:handles.Data1.nImages
            if min(handles.Data1.imageStack(i).data(:)) < Mn_o
                Mn_o = min(handles.Data1.imageStack(i).data(:));
            end
            if max(handles.Data1.imageStack(i).data(:)) > Mx_o
                Mx_o = max(handles.Data1.imageStack(i).data(:));
            end
        end
        
        handles.Current.clims1 = [Mn_o Mx_o];
        
        handles.Data1.clims_orig = handles.Current.clims1;
        set(handles.BlackSlider_C1,'Min',0,'Max',2*double(Mx_o),'SliderStep',[1/((2*double(Mx_o))-1),10/((2*double(Mx_o))-1)],'Value',double(Mn_o));
        set(handles.WhiteSlider_C1,'Min',0,'Max',2*double(Mx_o),'SliderStep',[1/((2*double(Mx_o))-1),10/((2*double(Mx_o))-1)],'Value',double(Mx_o));
        set(handles.BlackEdit_C1,'String',num2str(Mn_o));
        set(handles.WhiteEdit_C1,'String',num2str(Mx_o));
        if isfield(handles.Process1,'filterStack')
            set(handles.ImageViewPop,'Enable','on');
            
            Mn_f = 100000;
            Mx_f = 0;
            for i = 1:handles.Data1.nImages
                if min(handles.Process1.filterStack(i).data(:)) < Mn_f
                    Mn_f = min(handles.Process1.filterStack(i).data(:));
                end
                if max(handles.Process1.filterStack(i).data(:)) > Mx_f
                    Mx_f = max(handles.Process1.filterStack(i).data(:));
                end
                %             if Mx_f/Mn_f > 6
                %                 handles.Data1.clims_filt = [Mn_f 0.25*Mx_f];
                %             elseif Mx_o/Mn_o > 4
                %                 handles.Data1.clims_filt = [Mn_f 0.5*Mx_f];
                %             else
                handles.Data1.clims_filt = [Mn_f Mx_f];
                %             end
                
            end
        end
        %calculate the minimum & maximum values in the stack 2
        Mn_o = 100000;
        Mx_o = 0;
        
        for i = 1:handles.Data2.nImages
            if min(handles.Data2.imageStack(i).data(:)) < Mn_o
                Mn_o = min(handles.Data2.imageStack(i).data(:));
            end
            if max(handles.Data2.imageStack(i).data(:)) > Mx_o
                Mx_o = max(handles.Data2.imageStack(i).data(:));
            end
        end
        %
        handles.Current.clims2 = [Mn_o Mx_o];
        handles.Data2.clims_orig = handles.Current.clims2;
        set(handles.BlackSlider_C2,'Min',0,'Max',2*double(Mx_o),'SliderStep',[1/((2*double(Mx_o))-1),10/((2*double(Mx_o))-1)],'Value',double(Mn_o));
        set(handles.WhiteSlider_C2,'Min',0,'Max',2*double(Mx_o),'SliderStep',[1/((2*double(Mx_o))-1),10/((2*double(Mx_o))-1)],'Value',double(Mx_o));
        set(handles.BlackEdit_C2,'String',num2str(Mn_o));
        set(handles.WhiteEdit_C2,'String',num2str(Mx_o));
        
        %Enable other fields & plot if the data is available
        if isfield(handles.Process2,'filterStack')
            
            Mn_f = 100000;
            Mx_f = 0;
            for i = 1:handles.Data2.nImages
                if min(handles.Process2.filterStack(i).data(:)) < Mn_f
                    Mn_f = min(handles.Process2.filterStack(i).data(:));
                end
                if max(handles.Process2.filterStack(i).data(:)) > Mx_f
                    Mx_f = max(handles.Process2.filterStack(i).data(:));
                end
                %             if Mx_f/Mn_f > 6
                %                 handles.Data2.clims_filt = [Mn_f 0.25*Mx_f];
                %             elseif Mx_o/Mn_o > 4
                %                 handles.Data2.clims_filt = [Mn_f 0.5*Mx_f];
                %             else
                handles.Data2.clims_filt = [Mn_f Mx_f];
                %             end
                
            end
            
        end
        
        handles.Current.Image1 = ...
            handles.Data1.imageStack(handles.Current.Ind).data;
        
        axes(handles.Channel1data);
        colormap(gray);
        imagesc(handles.Current.Image1, handles.Current.clims1);
        axis image;
        set(handles.imageSlider,'Enable','on');
        % Set Minimum of image slider
        set(handles.imageSlider,'Min', 1);
        % Set Maximum of image slider
        set(handles.imageSlider,'Max', handles.Data1.nImages);
        % Set Value of image slider
        set(handles.imageSlider,'Value', handles.Current.Ind);
        % Set Step of image slider
        set(handles.imageSlider, 'SliderStep', [1/(handles.Data1.nImages-1)...
            1/(handles.Data1.nImages-1)]);
        % Set Image counter
        set(handles.ImageCount, 'String', ['Image ', ...
            num2str(handles.Current.Ind),'/',num2str(handles.Data1.nImages)]);
        
        if isfield(handles,'Data2')
            handles.Current.Image2 = handles.Data2.imageStack(handles.Current.Ind).data;
            axes(handles.Channel2data);
            colormap(gray);
            imagesc(handles.Current.Image2, handles.Current.clims2);
            axis image;
            
            %display the merged image
            %             if get(handles.MergeViewPop,'Value') == 1
            %                 if get(handles.MergeChannel1Check,'Value')
            %                     mergeChannel1 = handles.Current.Image1;
            %                 else
            %                     mergeChannel1 = zeros(size(handles.Current.Image1));
            %                 end
            %                 if get(handles.MergeChannel2Check,'Value')
            %                     mergeChannel2 = handles.Current.Image2;
            %                 else
            %                     mergeChannel2 = zeros(size(handles.Current.Image2));
            %                 end
            %                 handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
            handles.Current.mergedImage = merge2colorData(handles.Current.Image1,handles.Current.Image2,[1 2],handles.Current.clims1,handles.Current.clims2);
            handles.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,1),size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages);
            handles.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,1),size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages);
            for i = 1:handles.Data1.nImages
                handles.mergedStack.Image1(:,:,i) = handles.Data1.imageStack(i).data;
                handles.mergedStack.Image2(:,:,i) = handles.Data2.imageStack(i).data;
            end
            handles.Current.mergedStack = handles.mergedStack;
            handles.Current.MergedInd = handles.Current.Ind;
            set(handles.MergedData,'Visible','on');
            axes(handles.MergedData);
            image(handles.Current.mergedImage);
            set(handles.MergeChannel1Check,'Visible','on');
            set(handles.MergeChannel2Check,'Visible','on');
            set(handles.MergeChannel1Check,'Value',1);
            set(handles.MergeChannel2Check,'Value',1);
            set(handles.ViewSelectStatic,'Visible','on');
            set(handles.MergeViewPop,'Visible','on');
            set(handles.MergeViewPop,'Enable','on');
            set(handles.MergeViewPop,'Value',1);
            set(handles.ExportMergeImage,'Enable','on');
            axes(handles.MergedData);
            image(handles.Current.mergedImage);
            axis image
            %             end
        end
        if get(handles.showParticlesCheck, 'Value')
            if handles.isFitPSF1
                axes(handles.Channel1data);
                plotParticle(handles.Tracking1.Particles, ...
                    handles.Current.Ind, handles.isFitPSF1);
                
            else
                axes(handles.Channel1data);
                plotParticle(handles.Tracking1.Centroids, ...
                    handles.Current.Ind, handles.isFitPSF1);
            end
            
            if isfield(handles,'Tracking2') && isfield(handles,'isFitPSF2')
                if handles.isFitPSF2
                    axes(handles.Channel2data);
                    plotParticle(handles.Tracking2.Particles, ...
                        handles.Current.Ind, handles.isFitPSF2);
                    
                else
                    axes(handles.Channel2data);
                    plotParticle(handles.Tracking2.Centroids, ...
                        handles.Current.Ind, handles.isFitPSF2);
                end
            end
        end
        
        %Plot tracks if requested & available
        if get(handles.showTracksCheck, 'Value')
            axes(handles.Channel1data);
            plotTracks(handles.Current.Tracks1, handles.Current.Ind);
            if isfield(handles.Current,'Tracks2')
                axes(handles.Channel2data);
                plotTracks(handles.Current.Tracks2, handles.Current.Ind);
            end
        end
        
        if get(handles.showROIsCheck,'Value')
            if isfield(handles.Process1,'ROIpos')
                axes(handles.Channel1data);
                plotROI(handles.Process1.ROIpos);
            end
            if isfield(handles,'Process2');
                if isfield(handles.Process2,'ROIpos');
                    axes(handles.Channel2data);
                    plotROI(handles.Process2.ROIpos);
                end
            end
            
            
        end
        
        % Update the image indicator in the panel
        set(handles.ImageCount, 'String', ['Image ', num2str(handles.Current.Ind),'/',num2str(handles.Data1.nImages)]);
    end
    
end
guidata(hObject,handles);

% --- Executes on button press in ChanAlignPush.
function ChanAlignPush_Callback(hObject, eventdata, handles)
% hObject    handle to ChanAlignPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Images.Image1 = handles.Data1.imageStack;

Images.Image2 = handles.Data2.imageStack;

if handles.isFitPSF1 == 1
    Particles.Particles1 = handles.Tracking1.Particles;
else
    Particles.Particles1 = handles.Tracking1.Centroids;
end
if isfield(handles.Tracking1,'CheckTracks')
    
    Tracks.Tracks1 = handles.Tracking1.CheckTracks;
else
    
    Tracks.Tracks1 = handles.Tracking1.Tracks;
end

if handles.isFitPSF2 == 1
    Particles.Particles2 = handles.Tracking2.Particles;
else
    Particles.Particles2 = handles.Tracking2.Centroids;
end
if isfield(handles.Tracking2,'CheckTracks')
    
    Tracks.Tracks2 = handles.Tracking2.CheckTracks;
else
    
    Tracks.Tracks2 = handles.Tracking2.Tracks;
end
if handles.isAlign == 1
    errordlg('Alignment already performed. Process canceled','Alignment already exists');
else
    hChannelAlign = ChannelAlign_GUI(Images,Particles,Tracks);
    waitfor(hChannelAlign);
end





ImagesAlign.Image1 = getappdata(0,'imageStack1');
ImagesAlign.Image2 = getappdata(0,'imageStack2');
ParticlesAlign.Particles1 = getappdata(0,'Particles1');
ParticlesAlign.Particles2 = getappdata(0,'Particles2');
TracksAlign.Tracks1 = getappdata(0,'Tracks1');
TracksAlign.Tracks2 = getappdata(0,'Tracks2');



handles.isFitPSF1 = 1;
handles.isFitPSF2 = 1;

handles.isAlign = 1;
handles.Data1.imageStack = ImagesAlign.Image1;
handles.Data2.imageStack = ImagesAlign.Image2;
if handles.isFitPSF1 == 1
    handles.Tracking1.Particles = ParticlesAlign.Particles1;
else
    handles.Tracking1.Centroids = ParticlesAlign.Particles1;
end

if isfield(handles.Tracking1,'CheckTracks')
    
    handles.Tracking1.CheckTracks = TracksAlign.Tracks1;
else
    
    handles.Tracking1.Tracks = TracksAlign.Tracks1;
end

if handles.isFitPSF2 == 1
    handles.Tracking2.Particles = ParticlesAlign.Particles2;
else
    handles.Tracking2.Centroids = ParticlesAlign.Particles2;
end
if isfield(handles.Tracking2,'CheckTracks')
    
    handles.Tracking2.CheckTracks = TracksAlign.Tracks2;
else
    
    handles.Tracking2.Tracks = TracksAlign.Tracks2;
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in ColocPush.
function ColocPush_Callback(hObject, eventdata, handles)
% hObject    handle to ColocPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get parameters from user
defPar(1) = handles.Parameters.Bound.ThreshL;

defPar(2) = handles.Parameters.Bound.ThreshH;
defPar(3) = handles.Parameters.Bound.minBoundFrames;
defPar(4) = handles.Parameters.Used1.Acquisition.frameTime;

defaults = {num2str(defPar(1)), ... %num2str(defPar(2)), ...
    num2str(defPar(3)), '1', num2str(defPar(4))};
pxSize = str2double(get(handles.PxSizeText,'String'));
prompt = {'Maximum Jump between consecutive frames',...
    'Minimum length of bound tracks',...
    'Histogram binning(1 = no binning)',...
    'Frame Time'};
dlgtitle = 'Parameters for the residence time histogram';
AnalysisParam = inputdlg(prompt,dlgtitle, 1, defaults);
if ~isempty(AnalysisParam)
    handles.Parameters.Bound.ThreshL = str2double(AnalysisParam{1});
    % handles.Parameters.ThreshH = str2num(AnalysisS{2});
    handles.Parameters.Bound.minBoundFrames = str2double(AnalysisParam{2});
    handles.Parameters.bin = str2double(AnalysisParam{3});
    handles.Parameters.Used1.Acquisition.frameTime = str2double(AnalysisParam{4});
    handles.Parameters.Used2.Acquisition.frameTime = str2double(AnalysisParam{4});
    
    
    %Extract the total number of particles over time in order to fit the
    %photobleaching
    if isfield(handles.Tracking1,'Particles');
        Particles1 = handles.Tracking1.Particles;
        Particles1(:,10) = (Particles1(:,10) - 1)*pxSize;
        Particles1(:,11) = (Particles1(:,11) - 1)*pxSize;
    else
        Particles1 = handles.Tracking1.Centroids;
    end
    Particles1(:,1) = (Particles1(:,1) - 1)*pxSize;
    Particles1(:,2) = (Particles1(:,2) - 1)*pxSize;
    if isfield(handles.Tracking2,'Particles');
        Particles2 = handles.Tracking2.Particles;
        Particles2(:,10) = (Particles2(:,10) - 1)*pxSize;
        Particles2(:,11) = (Particles2(:,11) - 1)*pxSize;
    else
        Particles2 = handles.Tracking2.Centroids;
    end
    Particles2(:,1) = (Particles2(:,1) - 1)*pxSize;
    Particles2(:,2) = (Particles2(:,2) - 1)*pxSize;
    if isfield(handles.Tracking1,'CheckTracks') && isfield(handles.Tracking2,'CheckTracks')
        autoMan = questdlg('Do you want to analyze manually corrected tracks?','Manual or automatic?','Manual','Auto','Manual');
        if strcmp(autoMan,'Manual')
            T1 = handles.Tracking1.CheckTracks;
            T2 = handles.Tracking2.CheckTracks;
        else
            T1 = handles.Tracking1.Tracks;
            T2 = handles.Tracking2.Tracks;
        end
    else
        T1 = handles.Tracking1.Tracks;
        T2 = handles.Tracking2.Tracks;
    end
    
    
    T1(:,1:2) =	(T1(:,1:2) - 1).*pxSize;
    T2(:,1:2) =	(T2(:,1:2) - 1).*pxSize;
    
    
    ColocAnalysis(AnalysisParam,Particles1,Particles2,T1,T2);
    
end
% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function ExportMergeImage_Callback(hObject, eventdata, handles)
% hObject    handle to ExportMergeImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


MergeView = get(handles.MergeViewPop,'Value');

fileStr = 'Merged_Image';
if MergeView == 1;
    fileStr = [fileStr,'XY'];
elseif MergeView == 2
    fileStr = [fileStr,'XT'];
else
    fileStr = [fileStr,'YT'];
end

fileStr = [fileStr,'_',num2str(handles.Current.MergedInd)];

if ~isempty(handles.Data1.pathName)
    defName = [handles.Data1.pathName,fileStr];
else
    defName = [pwd,fileStr];
end

filtSpec = {'*.tif';'*.jpg'};

[fname, pname, filtIndex] = uiputfile(filtSpec,'Save Merged Image to File',defName);

if fname ~=0
    if filtIndex == 1
        fmt = 'tif';
    else
        fmt = 'jpeg';
    end
    imwrite(handles.Current.mergedImage,[pname,fname],fmt);
end




function FrameIntervalText_Callback(hObject, eventdata, handles)
% hObject    handle to FrameIntervalText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrameIntervalText as text
%        str2double(get(hObject,'String')) returns contents of FrameIntervalText as a double
Inter = str2double(get(hObject,'String'));
handles.Parameters.Acquisition1.FrameTime = Inter;
handles.Parameters.Acquisition2.FrameTime = Inter;

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function FrameIntervalText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameIntervalText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MergeColoc_Callback(hObject, eventdata, handles)
% hObject    handle to MergeColoc (see GCBO)
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
    FileNames = {FileNames};
end
% Set the parameters for the analysis of the residence time;
if ~isempty(FileNames)
    %get parameters from user
    defPar(1) = handles.Parameters.Bound.ThreshL;
    
    defPar(2) = handles.Parameters.Bound.ThreshH;
    defPar(3) = handles.Parameters.Bound.minBoundFrames;
    if isfield(handles.Parameters,'Used1')
        defPar(4) = handles.Parameters.Used1.Acquisition.frameTime;
    else
        defPar(4) = 0.2;
    end
    
    defaults = {num2str(defPar(1)), ... %num2str(defPar(2)), ...
        num2str(defPar(3)), '1', num2str(defPar(4)), '0'};
    pxSize = str2double(get(handles.PxSizeText,'String'));
    prompt = {'Maximum Jump between consecutive frames',...
        'Minimum length of bound tracks',...
        'Histogram binning(1 = no binning)',...
        'Frame Time',...
        'Max. Frames to Analyze (0 = all)'};
    dlgtitle = 'Parameters for the residence time histogram';
    AnalysisParam = inputdlg(prompt,dlgtitle, 1, defaults);
    if ~isempty(AnalysisParam)
        handles.Parameters.Bound.ThreshL = str2double(AnalysisParam{1});
        % handles.Parameters.ThreshH = str2num(AnalysisS{2});
        handles.Parameters.Bound.minBoundFrames = str2double(AnalysisParam{2});
        handles.Parameters.bin = str2double(AnalysisParam{3});
        handles.Parameters.Used1.Acquisition.frameTime = str2double(AnalysisParam{4});
        handles.Parameters.Used2.Acquisition.frameTime = str2double(AnalysisParam{4});
        handles.Parameters.maxAnalFrames = str2double(AnalysisParam{5});
        Particles1_all = [];
        Particles2_all = [];
        T1_all = [];
        T2_all = [];
        for i = 1:size(FileNames,2)
            if isempty(T1_all)
                trackNumStart1 = 0;
            else
                trackNumStart1 = max(T1_all(:,4));
            end
            if isempty(T2_all)
                trackNumStart2 = 0;
            else
                trackNumStart2 = max(T2_all(:,4));
            end
            Temp = load([PathName,FileNames{i}],'Results');
            offset = (i-1)*2000;
            %Extract the total number of particles over time in order to fit the
            %photobleaching
            if isfield(Temp.Results.Tracking1,'Particles')
                Particles1 = Temp.Results.Tracking1.Particles;
                if handles.Parameters.maxAnalFrames > 0
                    Particles1 = Particles1(Particles1(:,6) <= handles.Parameters.maxAnalFrames,:);
                end
                Particles1(:,10) = Particles1(:,10) + offset;
                Particles1(:,10) = (Particles1(:,10) - 1)*pxSize;
                Particles1(:,11) = (Particles1(:,11) - 1)*pxSize;
            else
                Particles1 = Temp.Results.Tracking1.Centroids;
            end
            Particles1(:,1) = Particles1(:,1) + offset;
            Particles1(:,1) = (Particles1(:,1) - 1)*pxSize;
            Particles1(:,2) = (Particles1(:,2) - 1)*pxSize;
            if isfield(Temp.Results.Tracking2,'Particles')
                Particles2 = Temp.Results.Tracking2.Particles;
                if handles.Parameters.maxAnalFrames > 0
                    Particles2 = Particles2(Particles2(:,6) <= handles.Parameters.maxAnalFrames,:);
                end
                Particles2(:,10) = Particles2(:,10) + offset;
                Particles2(:,10) = (Particles2(:,10) - 1)*pxSize;
                Particles2(:,11) = (Particles2(:,11) - 1)*pxSize;
            else
                Particles2 = Temp.Results.Tracking2.Centroids;
            end
            Particles2(:,1) = Particles2(:,1) + offset;
            Particles2(:,1) = (Particles2(:,1) - 1)*pxSize;
            Particles2(:,2) = (Particles2(:,2) - 1)*pxSize;
            if isfield(Temp.Results.Tracking1,'CheckTracks') && isfield(Temp.Results.Tracking2,'CheckTracks') && i == 1
                autoMan = questdlg('Do you want to analyze manually corrected tracks?','Manual or automatic?','Manual','Auto','Manual');
            elseif i == 1
                autoMan = 'Auto';
            end
            if strcmp(autoMan,'Manual')
                T1 = Temp.Results.Tracking1.CheckTracks;
                
                T2 = Temp.Results.Tracking2.CheckTracks;
            else
                T1 = Temp.Results.Tracking1.Tracks;
                T2 = Temp.Results.Tracking2.Tracks;
            end
            if handles.Parameters.maxAnalFrames > 0
                T1 = T1(T1(:,3) <= handles.Parameters.maxAnalFrames,:);
                T2 = T2(T2(:,3) <= handles.Parameters.maxAnalFrames,:);
            end
            T1(:,1) = T1(:,1) + offset;
            T1(:,4) = T1(:,4) + trackNumStart1;
            T2(:,1) = T2(:,1) + offset;
            T2(:,4) = T2(:,4) + trackNumStart2;
            
            T1(:,1:2) =	(T1(:,1:2) - 1).*pxSize;
            T2(:,1:2) =	(T2(:,1:2) - 1).*pxSize;
            
            
            Particles1_all = [Particles1_all; Particles1];
            Particles2_all = [Particles2_all; Particles2];
            T1_all = [T1_all; T1];
            T2_all = [T2_all; T2];
            
        end
        ColocAnalysis(AnalysisParam,Particles1_all,Particles2_all,T1_all,T2_all);
    end
    
end



function PxSizeText_Callback(hObject, eventdata, handles)
% hObject    handle to PxSizeText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PxSizeText as text
%        str2double(get(hObject,'String')) returns contents of PxSizeText as a double
handles.Parameters.Acquisition1.pixelSize = str2double(get(hObject,'String'));
handles.Parameters.Acquisition2.pixelSize = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function PxSizeText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PxSizeText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Auto_C1.
function Auto_C1_Callback(hObject, eventdata, handles)
% hObject    handle to Auto_C1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
what2draw = get(handles.ImageViewPop,'Value');

if what2draw == 1
    handles.Current.clims1 =  handles.Data1.clims_orig;
else
    handles.Current.clims1 =  handles.Data1.clims_filt;
end

set(handles.BlackEdit_C1,'String',num2str(handles.Current.clims1(1)));
set(handles.WhiteEdit_C1,'String',num2str(handles.Current.clims1(2)));

set(handles.BlackSlider_C1,'Val',handles.Current.clims1(1));
set(handles.WhiteSlider_C1,'Val',handles.Current.clims1(2));
set(handles.Channel1data,'NextPlot','replacechildren');
axes(handles.Channel1data);
imagesc(handles.Current.Image1, handles.Current.clims1);       % plot image
% axis image;
%Plot ROI

if get(handles.showParticlesCheck, 'Value')
    if handles.isFitPSF1
        plotParticle(handles.Tracking1.Particles, ...
            handles.Current.Ind, handles.isFitPSF1);
    else
        plotParticle(handles.Tracking1.Centroids, ...
            handles.Current.Ind, handles.isFitPSF1);
    end
end

if get(handles.showTracksCheck, 'Value')
    plotTracks(handles.Current.Tracks1, handles.Current.Ind);
end

if isfield(handles.Process1,'ROIpos') && get(handles.showROIsCheck,'Value')
    plotROI(handles.Process1.ROIpos);
    
end
%Show merged data

if get(handles.MergeViewPop,'Value') == 1
    if get(handles.MergeChannel1Check,'Value')
        mergeChannel1 = handles.Current.Image1;
    else
        mergeChannel1 = zeros(size(handles.Current.Image1));
    end
    if get(handles.MergeChannel2Check,'Value')
        mergeChannel2 = handles.Current.Image2;
    else
        mergeChannel2 = zeros(size(handles.Current.Image2));
    end
    handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    axis image
    
elseif get(handles.MergeViewPop,'Value') == 2 %X*T
    %         set(handles.mergeSlider,'Visible','on');
    %         set(handles.mergeSlider,'Enable','on');
    %
    %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,1));
    %
    %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,1));
    %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,1)-1) 1/(size(handles.Data1.imageStack(1).data,1)-1)]);
    %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[2 3 1]);
    %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[2 3 1]);
    %
    %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,1));
    %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,1));
    %         %         for i = 1:size(handles.Data2.imageStack(1).data,1)
    %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(i,:,:),size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages);
    %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(i,:,:),size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages);
    %         %         end
    %         handles.Current.MergedInd = 1;
    if get(handles.MergeChannel1Check,'Value')
        mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
    else
        mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
    end
    if get(handles.MergeChannel2Check,'Value')
        mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
    else
        mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
    end
    handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    xlabel('Frame Number');
    %         axis image
    set(handles.SliceCount,'Visible','on');
    set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,1))]);
    set(handles.LocationDrawCheck,'Visible','on');
    %         ax(3) = [];
    %         linkaxes(ax);
    
elseif get(handles.MergeViewPop,'Value') == 3 %Y*T
    %         set(handles.mergeSlider,'Visible','on');
    %         set(handles.mergeSlider,'Enable','on');
    %
    %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,2))
    %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,2));
    %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,2)-1) 1/(size(handles.Data1.imageStack(1).data,2)-1)]);
    %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[1 3 2]);
    %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[1 3 2]);
    %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,2));
    %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,2));
    %         %         for i = 1:size(handles.Data2.imageStack(1).data,2)
    %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(:,i,:),size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages);
    %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(:,i,:),size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages);
    %         %         end
    %         handles.Current.MergedInd = 1;
    if get(handles.MergeChannel1Check,'Value')
        mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
    else
        mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
    end
    if get(handles.MergeChannel2Check,'Value')
        mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
    else
        mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
    end
    handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    xlabel('Frame Number');
    %         axis image
    set(handles.SliceCount,'Visible','on');
    set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,2))]);
    set(handles.LocationDrawCheck,'Visible','on');
    %         ax(3) = [];
    %         linkaxes(ax);
end

if get(handles.LocationDrawCheck,'Value')
    if get(handles.MergeViewPop,'Value') == 2
        tyline = [0 size(handles.Current.Image1,2)];
        txline = [handles.Current.Ind handles.Current.Ind];
        xline = [0 size(handles.Current.Image1,2)];
        yline = [handles.Current.MergedInd handles.Current.MergedInd];
    elseif get(handles.MergeViewPop,'Value') == 3
        tyline = [0 size(handles.Current.Image1,1)];
        txline = [handles.Current.Ind handles.Current.Ind];
        yline = [0 size(handles.Current.Image1,1)];
        xline = [handles.Current.MergedInd handles.Current.MergedInd];
    end
    axes(handles.Channel1data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    axes(handles.Channel2data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    
    axes(handles.MergedData);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(txline,tyline,'w','LineWidth',0.75);
    hold off
end
guidata(hObject,handles);



% --- Executes on slider movement.
function WhiteSlider_C1_Callback(hObject, eventdata, handles)
% hObject    handle to WhiteSlider_C1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
MaxBrightnessStack = get(handles.BlackSlider_C1,'Max');

WhiteLevel_C1 = round(double((get(hObject,'Value')))*100)/100;
%Make sure we don't try to set the brightness outside of the slider limits
if WhiteLevel_C1 <= 0
    WhiteLevel_C1 = 0.1;
elseif WhiteLevel_C1 > MaxBrightnessStack
    WhiteLevel_C1 = MaxBrightnessStack-0.1;
end
set(handles.WhiteEdit_C1,'String',num2str(WhiteLevel_C1));

BlackLevel_C1 = str2double(get(handles.BlackEdit_C1,'String'));

if BlackLevel_C1 >= WhiteLevel_C1
    BlackLevel_C1 = WhiteLevel_C1-0.1;
    set(handles.BlackEdit_C1,'String',num2str(round(double(BlackLevel_C1)*100)/100));
end

set(handles.BlackSlider_C1,'Val',BlackLevel_C1);
set(handles.WhiteSlider_C1,'Val',WhiteLevel_C1);
handles.Current.clims1 = [BlackLevel_C1 WhiteLevel_C1];


set(handles.Channel1data,'NextPlot','replacechildren');
axes(handles.Channel1data);
imagesc(handles.Current.Image1, handles.Current.clims1);       % plot image
% axis image;
%Plot ROI

if get(handles.showParticlesCheck, 'Value')
    if handles.isFitPSF1
        plotParticle(handles.Tracking1.Particles, ...
            handles.Current.Ind, handles.isFitPSF1);
    else
        plotParticle(handles.Tracking1.Centroids, ...
            handles.Current.Ind, handles.isFitPSF1);
    end
end

if get(handles.showTracksCheck, 'Value')
    plotTracks(handles.Current.Tracks1, handles.Current.Ind);
end

if isfield(handles.Process1,'ROIpos') && get(handles.showROIsCheck,'Value')
    plotROI(handles.Process1.ROIpos);
    
end
%Show merged data
% if isfield(handles,'MergedInd')

    if get(handles.MergeViewPop,'Value') == 1
        if get(handles.MergeChannel1Check,'Value')
            mergeChannel1 = handles.Current.Image1;
        else
            mergeChannel1 = zeros(size(handles.Current.Image1));
        end
        if get(handles.MergeChannel2Check,'Value')
            mergeChannel2 = handles.Current.Image2;
        else
            mergeChannel2 = zeros(size(handles.Current.Image2));
        end
        handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
        axes(handles.MergedData);
        image(handles.Current.mergedImage);
        axis image

    elseif get(handles.MergeViewPop,'Value') == 2 %X*T
        %         set(handles.mergeSlider,'Visible','on');
        %         set(handles.mergeSlider,'Enable','on');
        %
        %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,1));
        %
        %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,1));
        %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,1)-1) 1/(size(handles.Data1.imageStack(1).data,1)-1)]);
        %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[2 3 1]);
        %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[2 3 1]);
        %
        %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,1));
        %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,1));
        %         %         for i = 1:size(handles.Data2.imageStack(1).data,1)
        %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(i,:,:),size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages);
        %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(i,:,:),size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages);
        %         %         end
        %         handles.Current.MergedInd = 1;
        if get(handles.MergeChannel1Check,'Value')
            mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
        else
            mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
        end
        if get(handles.MergeChannel2Check,'Value')
            mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
        else
            mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
        end
        handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
        axes(handles.MergedData);
        image(handles.Current.mergedImage);
        xlabel('Frame Number');
        %         axis image
        set(handles.SliceCount,'Visible','on');
        set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,1))]);
        set(handles.LocationDrawCheck,'Visible','on');
        %         ax(3) = [];
        %         linkaxes(ax);

    elseif get(handles.MergeViewPop,'Value') == 3 %Y*T
        %         set(handles.mergeSlider,'Visible','on');
        %         set(handles.mergeSlider,'Enable','on');
        %
        %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,2))
        %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,2));
        %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,2)-1) 1/(size(handles.Data1.imageStack(1).data,2)-1)]);
        %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[1 3 2]);
        %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[1 3 2]);
        %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,2));
        %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,2));
        %         %         for i = 1:size(handles.Data2.imageStack(1).data,2)
        %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(:,i,:),size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages);
        %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(:,i,:),size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages);
        %         %         end
        %         handles.Current.MergedInd = 1;
        if get(handles.MergeChannel1Check,'Value')
            mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
        else
            mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
        end
        if get(handles.MergeChannel2Check,'Value')
            mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
        else
            mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
        end
        handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
        axes(handles.MergedData);
        image(handles.Current.mergedImage);
        xlabel('Frame Number');
        %         axis image
        set(handles.SliceCount,'Visible','on');
        set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,2))]);
        set(handles.LocationDrawCheck,'Visible','on');
        %         ax(3) = [];
        %         linkaxes(ax);
    end
% end
if get(handles.LocationDrawCheck,'Value')
    if get(handles.MergeViewPop,'Value') == 2
        tyline = [0 size(handles.Current.Image1,2)];
        txline = [handles.Current.Ind handles.Current.Ind];
        xline = [0 size(handles.Current.Image1,2)];
        yline = [handles.Current.MergedInd handles.Current.MergedInd];
    elseif get(handles.MergeViewPop,'Value') == 3
        tyline = [0 size(handles.Current.Image1,1)];
        txline = [handles.Current.Ind handles.Current.Ind];
        yline = [0 size(handles.Current.Image1,1)];
        xline = [handles.Current.MergedInd handles.Current.MergedInd];
    end
    axes(handles.Channel1data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    axes(handles.Channel2data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    
    axes(handles.MergedData);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(txline,tyline,'w','LineWidth',0.75);
    hold off
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function WhiteSlider_C1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WhiteSlider_C1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function BlackSlider_C1_Callback(hObject, eventdata, handles)
% hObject    handle to BlackSlider_C1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
MaxBrightnessStack = get(handles.BlackSlider_C1,'Max');

BlackLevel_C1 = round(double((get(hObject,'Value')))*100)/100;
%Make sure we don't try to set the brightness outside of the slider limits
if BlackLevel_C1 < 0
    BlackLevel_C1 = 0;
elseif BlackLevel_C1 > MaxBrightnessStack
    BlackLevel_C1 = MaxBrightnessStack-0.1;
end
set(handles.BlackEdit_C1,'String',num2str(BlackLevel_C1));

WhiteLevel_C1 = str2double(get(handles.WhiteEdit_C1,'String'));

if BlackLevel_C1 >= WhiteLevel_C1
    WhiteLevel_C1 = BlackLevel_C1+0.1;
    set(handles.WhiteEdit_C1,'String',num2str(round(double(WhiteLevel_C1)*100)/100));
end

set(handles.BlackSlider_C1,'Val',BlackLevel_C1);
set(handles.WhiteSlider_C1,'Val',WhiteLevel_C1);
handles.Current.clims1 = [BlackLevel_C1 WhiteLevel_C1];


set(handles.Channel1data,'NextPlot','replacechildren');
axes(handles.Channel1data);
imagesc(handles.Current.Image1, handles.Current.clims1);       % plot image
% axis image;
%Plot ROI

if get(handles.showParticlesCheck, 'Value')
    if handles.isFitPSF1
        plotParticle(handles.Tracking1.Particles, ...
            handles.Current.Ind, handles.isFitPSF1);
    else
        plotParticle(handles.Tracking1.Centroids, ...
            handles.Current.Ind, handles.isFitPSF1);
    end
end

if get(handles.showTracksCheck, 'Value')
    plotTracks(handles.Current.Tracks1, handles.Current.Ind);
end

if isfield(handles.Process1,'ROIpos') && get(handles.showROIsCheck,'Value')
    plotROI(handles.Process1.ROIpos);
    
end
%Show merged data
% if isfield(handles,'MergedInd')

    if get(handles.MergeViewPop,'Value') == 1
        if get(handles.MergeChannel1Check,'Value')
            mergeChannel1 = handles.Current.Image1;
        else
            mergeChannel1 = zeros(size(handles.Current.Image1));
        end
        if get(handles.MergeChannel2Check,'Value')
            mergeChannel2 = handles.Current.Image2;
        else
            mergeChannel2 = zeros(size(handles.Current.Image2));
        end
        handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
        axes(handles.MergedData);
        image(handles.Current.mergedImage);
        axis image

    elseif get(handles.MergeViewPop,'Value') == 2 %X*T
        %         set(handles.mergeSlider,'Visible','on');
        %         set(handles.mergeSlider,'Enable','on');
        %
        %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,1));
        %
        %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,1));
        %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,1)-1) 1/(size(handles.Data1.imageStack(1).data,1)-1)]);
        %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[2 3 1]);
        %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[2 3 1]);
        %
        %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,1));
        %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,1));
        %         %         for i = 1:size(handles.Data2.imageStack(1).data,1)
        %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(i,:,:),size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages);
        %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(i,:,:),size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages);
        %         %         end
        %         handles.Current.MergedInd = 1;
        if get(handles.MergeChannel1Check,'Value')
            mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
        else
            mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
        end
        if get(handles.MergeChannel2Check,'Value')
            mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
        else
            mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
        end
        handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
        axes(handles.MergedData);
        image(handles.Current.mergedImage);
        xlabel('Frame Number');
        %         axis image
        set(handles.SliceCount,'Visible','on');
        set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,1))]);
        set(handles.LocationDrawCheck,'Visible','on');
        %         ax(3) = [];
        %         linkaxes(ax);

    elseif get(handles.MergeViewPop,'Value') == 3 %Y*T
        %         set(handles.mergeSlider,'Visible','on');
        %         set(handles.mergeSlider,'Enable','on');
        %
        %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,2))
        %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,2));
        %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,2)-1) 1/(size(handles.Data1.imageStack(1).data,2)-1)]);
        %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[1 3 2]);
        %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[1 3 2]);
        %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,2));
        %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,2));
        %         %         for i = 1:size(handles.Data2.imageStack(1).data,2)
        %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(:,i,:),size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages);
        %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(:,i,:),size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages);
        %         %         end
        %         handles.Current.MergedInd = 1;
        if get(handles.MergeChannel1Check,'Value')
            mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
        else
            mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
        end
        if get(handles.MergeChannel2Check,'Value')
            mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
        else
            mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
        end
        handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
        axes(handles.MergedData);
        image(handles.Current.mergedImage);
        xlabel('Frame Number');
        %         axis image
        set(handles.SliceCount,'Visible','on');
        set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,2))]);
        set(handles.LocationDrawCheck,'Visible','on');
        %         ax(3) = [];
        %         linkaxes(ax);
    end
% end
if get(handles.LocationDrawCheck,'Value')
    if get(handles.MergeViewPop,'Value') == 2
        tyline = [0 size(handles.Current.Image1,2)];
        txline = [handles.Current.Ind handles.Current.Ind];
        xline = [0 size(handles.Current.Image1,2)];
        yline = [handles.Current.MergedInd handles.Current.MergedInd];
    elseif get(handles.MergeViewPop,'Value') == 3
        tyline = [0 size(handles.Current.Image1,1)];
        txline = [handles.Current.Ind handles.Current.Ind];
        yline = [0 size(handles.Current.Image1,1)];
        xline = [handles.Current.MergedInd handles.Current.MergedInd];
    end
    axes(handles.Channel1data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    axes(handles.Channel2data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    
    axes(handles.MergedData);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(txline,tyline,'w','LineWidth',0.75);
    hold off
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function BlackSlider_C1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BlackSlider_C1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function BlackEdit_C1_Callback(hObject, eventdata, handles)
% hObject    handle to BlackEdit_C1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BlackEdit_C1 as text
%        str2double(get(hObject,'String')) returns contents of BlackEdit_C1 as a double
MaxBrightnessStack = get(handles.BlackSlider_C1,'Max');

BlackLevel_C1 = str2double(get(hObject,'String'));
%Make sure we don't try to set the brightness outside of the slider limits
if BlackLevel_C1 < 0
    BlackLevel_C1 = 0;
elseif BlackLevel_C1 > MaxBrightnessStack
    BlackLevel_C1 = MaxBrightnessStack-0.1;
end

set(hObject,'String',num2str(round(double(BlackLevel_C1)*100)/100));
WhiteLevel_C1 = str2double(get(handles.WhiteEdit_C1,'String'));

if BlackLevel_C1 >= WhiteLevel_C1
    WhiteLevel_C1 = BlackLevel_C1+0.1;
    set(handles.WhiteEdit_C1,'String',num2str(round(double(WhiteLevel_C1)*100)/100));
end

set(handles.BlackSlider_C1,'Val',BlackLevel_C1);
set(handles.WhiteSlider_C1,'Val',WhiteLevel_C1);
handles.Current.clims1 = [BlackLevel_C1 WhiteLevel_C1];


set(handles.Channel1data,'NextPlot','replacechildren');
axes(handles.Channel1data);
imagesc(handles.Current.Image1, handles.Current.clims1);       % plot image
% axis image;
%Plot ROI

if get(handles.showParticlesCheck, 'Value')
    if handles.isFitPSF1
        plotParticle(handles.Tracking1.Particles, ...
            handles.Current.Ind, handles.isFitPSF1);
    else
        plotParticle(handles.Tracking1.Centroids, ...
            handles.Current.Ind, handles.isFitPSF1);
    end
end

if get(handles.showTracksCheck, 'Value')
    plotTracks(handles.Current.Tracks1, handles.Current.Ind);
end

if isfield(handles.Process1,'ROIpos') && get(handles.showROIsCheck,'Value')
    plotROI(handles.Process1.ROIpos);
    
end
%Show merged data

if get(handles.MergeViewPop,'Value') == 1
    if get(handles.MergeChannel1Check,'Value')
        mergeChannel1 = handles.Current.Image1;
    else
        mergeChannel1 = zeros(size(handles.Current.Image1));
    end
    if get(handles.MergeChannel2Check,'Value')
        mergeChannel2 = handles.Current.Image2;
    else
        mergeChannel2 = zeros(size(handles.Current.Image2));
    end
    handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    axis image
    
elseif get(handles.MergeViewPop,'Value') == 2 %X*T
    %         set(handles.mergeSlider,'Visible','on');
    %         set(handles.mergeSlider,'Enable','on');
    %
    %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,1));
    %
    %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,1));
    %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,1)-1) 1/(size(handles.Data1.imageStack(1).data,1)-1)]);
    %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[2 3 1]);
    %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[2 3 1]);
    %
    %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,1));
    %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,1));
    %         %         for i = 1:size(handles.Data2.imageStack(1).data,1)
    %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(i,:,:),size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages);
    %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(i,:,:),size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages);
    %         %         end
    %         handles.Current.MergedInd = 1;
    if get(handles.MergeChannel1Check,'Value')
        mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
    else
        mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
    end
    if get(handles.MergeChannel2Check,'Value')
        mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
    else
        mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
    end
    handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    xlabel('Frame Number');
    %         axis image
    set(handles.SliceCount,'Visible','on');
    set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,1))]);
    set(handles.LocationDrawCheck,'Visible','on');
    %         ax(3) = [];
    %         linkaxes(ax);
    
elseif get(handles.MergeViewPop,'Value') == 3 %Y*T
    %         set(handles.mergeSlider,'Visible','on');
    %         set(handles.mergeSlider,'Enable','on');
    %
    %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,2))
    %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,2));
    %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,2)-1) 1/(size(handles.Data1.imageStack(1).data,2)-1)]);
    %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[1 3 2]);
    %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[1 3 2]);
    %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,2));
    %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,2));
    %         %         for i = 1:size(handles.Data2.imageStack(1).data,2)
    %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(:,i,:),size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages);
    %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(:,i,:),size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages);
    %         %         end
    %         handles.Current.MergedInd = 1;
    if get(handles.MergeChannel1Check,'Value')
        mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
    else
        mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
    end
    if get(handles.MergeChannel2Check,'Value')
        mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
    else
        mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
    end
    handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    xlabel('Frame Number');
    %         axis image
    set(handles.SliceCount,'Visible','on');
    set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,2))]);
    set(handles.LocationDrawCheck,'Visible','on');
    %         ax(3) = [];
    %         linkaxes(ax);
end

if get(handles.LocationDrawCheck,'Value')
    if get(handles.MergeViewPop,'Value') == 2
        tyline = [0 size(handles.Current.Image1,2)];
        txline = [handles.Current.Ind handles.Current.Ind];
        xline = [0 size(handles.Current.Image1,2)];
        yline = [handles.Current.MergedInd handles.Current.MergedInd];
    elseif get(handles.MergeViewPop,'Value') == 3
        tyline = [0 size(handles.Current.Image1,1)];
        txline = [handles.Current.Ind handles.Current.Ind];
        yline = [0 size(handles.Current.Image1,1)];
        xline = [handles.Current.MergedInd handles.Current.MergedInd];
    end
    axes(handles.Channel1data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    axes(handles.Channel2data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    
    axes(handles.MergedData);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(txline,tyline,'w','LineWidth',0.75);
    hold off
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function BlackEdit_C1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BlackEdit_C1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function WhiteEdit_C1_Callback(hObject, eventdata, handles)
% hObject    handle to WhiteEdit_C1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WhiteEdit_C1 as text
%        str2double(get(hObject,'String')) returns contents of WhiteEdit_C1 as a double
MaxBrightnessStack = get(handles.BlackSlider_C1,'Max');

WhiteLevel_C1 = str2double(get(hObject,'String'));
%Make sure we don't try to set the brightness outside of the slider limits
if WhiteLevel_C1 <= 0
    WhiteLevel_C1 = 0.1;
elseif WhiteLevel_C1 > MaxBrightnessStack
    WhiteLevel_C1 = MaxBrightnessStack;
end

set(hObject,'String',num2str(round(double(WhiteLevel_C1)*100)/100));
BlackLevel_C1 = str2double(get(handles.BlackEdit_C1,'String'));

if BlackLevel_C1 >= WhiteLevel_C1
    BlackLevel_C1 = WhiteLevel_C1-0.1;
    set(handles.BlackEdit_C1,'String',num2str(round(double(BlackLevel_C1)*100)/100));
end

set(handles.BlackSlider_C1,'Val',BlackLevel_C1);
set(handles.WhiteSlider_C1,'Val',WhiteLevel_C1);
handles.Current.clims1 = [BlackLevel_C1 WhiteLevel_C1];


set(handles.Channel1data,'NextPlot','replacechildren');
axes(handles.Channel1data);
imagesc(handles.Current.Image1, handles.Current.clims1);       % plot image
% axis image;
%Plot ROI

if get(handles.showParticlesCheck, 'Value')
    if handles.isFitPSF1
        plotParticle(handles.Tracking1.Particles, ...
            handles.Current.Ind, handles.isFitPSF1);
    else
        plotParticle(handles.Tracking1.Centroids, ...
            handles.Current.Ind, handles.isFitPSF1);
    end
end

if get(handles.showTracksCheck, 'Value')
    plotTracks(handles.Current.Tracks1, handles.Current.Ind);
end

if isfield(handles.Process1,'ROIpos') && get(handles.showROIsCheck,'Value')
    plotROI(handles.Process1.ROIpos);
    
end

if get(handles.MergeViewPop,'Value') == 1
    if get(handles.MergeChannel1Check,'Value')
        mergeChannel1 = handles.Current.Image1;
    else
        mergeChannel1 = zeros(size(handles.Current.Image1));
    end
    if get(handles.MergeChannel2Check,'Value')
        mergeChannel2 = handles.Current.Image2;
    else
        mergeChannel2 = zeros(size(handles.Current.Image2));
    end
    handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    axis image
elseif get(handles.MergeViewPop,'Value') == 2 %X*T
    %         set(handles.mergeSlider,'Visible','on');
    %         set(handles.mergeSlider,'Enable','on');
    %
    %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,1));
    %
    %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,1));
    %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,1)-1) 1/(size(handles.Data1.imageStack(1).data,1)-1)]);
    %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[2 3 1]);
    %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[2 3 1]);
    %
    %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,1));
    %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,1));
    %         %         for i = 1:size(handles.Data2.imageStack(1).data,1)
    %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(i,:,:),size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages);
    %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(i,:,:),size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages);
    %         %         end
    %         handles.Current.MergedInd = 1;
    if get(handles.MergeChannel1Check,'Value')
        mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
    else
        mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
    end
    if get(handles.MergeChannel2Check,'Value')
        mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
    else
        mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
    end
    handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    xlabel('Frame Number');
    %         axis image
    set(handles.SliceCount,'Visible','on');
    set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,1))]);
    set(handles.LocationDrawCheck,'Visible','on');
    %         ax(3) = [];
    %         linkaxes(ax);
    
elseif get(handles.MergeViewPop,'Value') == 3 %Y*T
    %         set(handles.mergeSlider,'Visible','on');
    %         set(handles.mergeSlider,'Enable','on');
    %
    %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,2))
    %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,2));
    %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,2)-1) 1/(size(handles.Data1.imageStack(1).data,2)-1)]);
    %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[1 3 2]);
    %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[1 3 2]);
    %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,2));
    %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,2));
    %         %         for i = 1:size(handles.Data2.imageStack(1).data,2)
    %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(:,i,:),size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages);
    %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(:,i,:),size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages);
    %         %         end
    %         handles.Current.MergedInd = 1;
    if get(handles.MergeChannel1Check,'Value')
        mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
    else
        mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
    end
    if get(handles.MergeChannel2Check,'Value')
        mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
    else
        mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
    end
    handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    xlabel('Frame Number');
    %         axis image
    set(handles.SliceCount,'Visible','on');
    set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,2))]);
    set(handles.LocationDrawCheck,'Visible','on');
    %         ax(3) = [];
    %         linkaxes(ax);
end
if get(handles.LocationDrawCheck,'Value')
    if get(handles.MergeViewPop,'Value') == 2
        tyline = [0 size(handles.Current.Image1,2)];
        txline = [handles.Current.Ind handles.Current.Ind];
        xline = [0 size(handles.Current.Image1,2)];
        yline = [handles.Current.MergedInd handles.Current.MergedInd];
    elseif get(handles.MergeViewPop,'Value') == 3
        tyline = [0 size(handles.Current.Image1,1)];
        txline = [handles.Current.Ind handles.Current.Ind];
        yline = [0 size(handles.Current.Image1,1)];
        xline = [handles.Current.MergedInd handles.Current.MergedInd];
    end
    axes(handles.Channel1data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    axes(handles.Channel2data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    
    axes(handles.MergedData);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(txline,tyline,'w','LineWidth',0.75);
    hold off
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function WhiteEdit_C1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WhiteEdit_C1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Auto_C2.
function Auto_C2_Callback(hObject, eventdata, handles)
% hObject    handle to Auto_C2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
what2draw = get(handles.ImageViewPop,'Value');

if what2draw == 1
    handles.Current.clims2 =  handles.Data2.clims_orig;
else
    handles.Current.clims2 =  handles.Data2.clims_filt;
end

set(handles.Channel2data,'NextPlot','replacechildren');
axes(handles.Channel2data);
imagesc(handles.Current.Image2, handles.Current.clims2);       % plot image
% axis image;
%Plot ROI

if get(handles.showParticlesCheck, 'Value')
    if handles.isFitPSF2
        plotParticle(handles.Tracking2.Particles, ...
            handles.Current.Ind, handles.isFitPSF2);
    else
        plotParticle(handles.Tracking2.Centroids, ...
            handles.Current.Ind, handles.isFitPSF2);
    end
end

if get(handles.showTracksCheck, 'Value')
    plotTracks(handles.Current.Tracks2, handles.Current.Ind);
end

if isfield(handles.Process2,'ROIpos') && get(handles.showROIsCheck,'Value')
    plotROI(handles.Process2.ROIpos);
    
end


if get(handles.MergeViewPop,'Value') == 1
    if get(handles.MergeChannel1Check,'Value')
        mergeChannel1 = handles.Current.Image1;
    else
        mergeChannel1 = zeros(size(handles.Current.Image1));
    end
    if get(handles.MergeChannel2Check,'Value')
        mergeChannel2 = handles.Current.Image2;
    else
        mergeChannel2 = zeros(size(handles.Current.Image2));
    end
    handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    axis image
elseif get(handles.MergeViewPop,'Value') == 2 %X*T
    %         set(handles.mergeSlider,'Visible','on');
    %         set(handles.mergeSlider,'Enable','on');
    %
    %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,1));
    %
    %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,1));
    %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,1)-1) 1/(size(handles.Data1.imageStack(1).data,1)-1)]);
    %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[2 3 1]);
    %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[2 3 1]);
    %
    %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,1));
    %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,1));
    %         %         for i = 1:size(handles.Data2.imageStack(1).data,1)
    %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(i,:,:),size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages);
    %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(i,:,:),size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages);
    %         %         end
    %         handles.Current.MergedInd = 1;
    if get(handles.MergeChannel1Check,'Value')
        mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
    else
        mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
    end
    if get(handles.MergeChannel2Check,'Value')
        mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
    else
        mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
    end
    handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    xlabel('Frame Number');
    %         axis image
    set(handles.SliceCount,'Visible','on');
    set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,1))]);
    set(handles.LocationDrawCheck,'Visible','on');
    %         ax(3) = [];
    %         linkaxes(ax);
    
elseif get(handles.MergeViewPop,'Value') == 3 %Y*T
    %         set(handles.mergeSlider,'Visible','on');
    %         set(handles.mergeSlider,'Enable','on');
    %
    %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,2))
    %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,2));
    %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,2)-1) 1/(size(handles.Data1.imageStack(1).data,2)-1)]);
    %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[1 3 2]);
    %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[1 3 2]);
    %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,2));
    %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,2));
    %         %         for i = 1:size(handles.Data2.imageStack(1).data,2)
    %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(:,i,:),size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages);
    %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(:,i,:),size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages);
    %         %         end
    %         handles.Current.MergedInd = 1;
    if get(handles.MergeChannel1Check,'Value')
        mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
    else
        mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
    end
    if get(handles.MergeChannel2Check,'Value')
        mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
    else
        mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
    end
    handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    xlabel('Frame Number');
    %         axis image
    set(handles.SliceCount,'Visible','on');
    set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,2))]);
    set(handles.LocationDrawCheck,'Visible','on');
    %         ax(3) = [];
    %         linkaxes(ax);
end

if get(handles.LocationDrawCheck,'Value')
    if get(handles.MergeViewPop,'Value') == 2
        tyline = [0 size(handles.Current.Image1,2)];
        txline = [handles.Current.Ind handles.Current.Ind];
        xline = [0 size(handles.Current.Image1,2)];
        yline = [handles.Current.MergedInd handles.Current.MergedInd];
    elseif get(handles.MergeViewPop,'Value') == 3
        tyline = [0 size(handles.Current.Image1,1)];
        txline = [handles.Current.Ind handles.Current.Ind];
        yline = [0 size(handles.Current.Image1,1)];
        xline = [handles.Current.MergedInd handles.Current.MergedInd];
    end
    axes(handles.Channel1data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    axes(handles.Channel2data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    
    axes(handles.MergedData);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(txline,tyline,'w','LineWidth',0.75);
    hold off
end
guidata(hObject,handles);

% --- Executes on slider movement.
function WhiteSlider_C2_Callback(hObject, eventdata, handles)
% hObject    handle to WhiteSlider_C2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
MaxBrightnessStack = get(handles.BlackSlider_C2,'Max');

WhiteLevel_C2 = round(double((get(hObject,'Value')))*100)/100;
%Make sure we don't try to set the brightness outside of the slider limits
if WhiteLevel_C2 <= 0
    WhiteLevel_C2 = 0.1;
elseif WhiteLevel_C2 > MaxBrightnessStack
    WhiteLevel_C2 = MaxBrightnessStack-0.1;
end
set(handles.WhiteEdit_C2,'String',num2str(WhiteLevel_C2));

BlackLevel_C2 = str2double(get(handles.BlackEdit_C2,'String'));

if BlackLevel_C2 >= WhiteLevel_C2
    BlackLevel_C2 = WhiteLevel_C2-0.1;
    set(handles.BlackEdit_C2,'String',num2str(round(double(BlackLevel_C2)*100)/100));
end

set(handles.BlackSlider_C2,'Val',BlackLevel_C2);
set(handles.WhiteSlider_C2,'Val',WhiteLevel_C2);
handles.Current.clims2 = [BlackLevel_C2 WhiteLevel_C2];


set(handles.Channel2data,'NextPlot','replacechildren');
axes(handles.Channel2data);
imagesc(handles.Current.Image2, handles.Current.clims2);       % plot image
% axis image;
%Plot ROI

if get(handles.showParticlesCheck, 'Value')
    if handles.isFitPSF2
        plotParticle(handles.Tracking2.Particles, ...
            handles.Current.Ind, handles.isFitPSF2);
    else
        plotParticle(handles.Tracking2.Centroids, ...
            handles.Current.Ind, handles.isFitPSF2);
    end
end

if get(handles.showTracksCheck, 'Value')
    plotTracks(handles.Current.Tracks2, handles.Current.Ind);
end

if isfield(handles.Process2,'ROIpos') && get(handles.showROIsCheck,'Value')
    plotROI(handles.Process2.ROIpos);
    
end
%Show merged data
% if isfield(handles,'MergedInd')

    if get(handles.MergeViewPop,'Value') == 1
        if get(handles.MergeChannel1Check,'Value')
            mergeChannel1 = handles.Current.Image1;
        else
            mergeChannel1 = zeros(size(handles.Current.Image1));
        end
        if get(handles.MergeChannel2Check,'Value')
            mergeChannel2 = handles.Current.Image2;
        else
            mergeChannel2 = zeros(size(handles.Current.Image2));
        end
        handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
        axes(handles.MergedData);
        image(handles.Current.mergedImage);
        axis image

    elseif get(handles.MergeViewPop,'Value') == 2 %X*T
        %         set(handles.mergeSlider,'Visible','on');
        %         set(handles.mergeSlider,'Enable','on');
        %
        %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,1));
        %
        %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,1));
        %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,1)-1) 1/(size(handles.Data1.imageStack(1).data,1)-1)]);
        %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[2 3 1]);
        %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[2 3 1]);
        %
        %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,1));
        %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,1));
        %         %         for i = 1:size(handles.Data2.imageStack(1).data,1)
        %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(i,:,:),size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages);
        %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(i,:,:),size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages);
        %         %         end
        %         handles.Current.MergedInd = 1;
        if get(handles.MergeChannel1Check,'Value')
            mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
        else
            mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
        end
        if get(handles.MergeChannel2Check,'Value')
            mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
        else
            mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
        end
        handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
        axes(handles.MergedData);
        image(handles.Current.mergedImage);
        xlabel('Frame Number');
        %         axis image
        set(handles.SliceCount,'Visible','on');
        set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,1))]);
        set(handles.LocationDrawCheck,'Visible','on');
        %         ax(3) = [];
        %         linkaxes(ax);

    elseif get(handles.MergeViewPop,'Value') == 3 %Y*T
        %         set(handles.mergeSlider,'Visible','on');
        %         set(handles.mergeSlider,'Enable','on');
        %
        %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,2))
        %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,2));
        %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,2)-1) 1/(size(handles.Data1.imageStack(1).data,2)-1)]);
        %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[1 3 2]);
        %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[1 3 2]);
        %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,2));
        %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,2));
        %         %         for i = 1:size(handles.Data2.imageStack(1).data,2)
        %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(:,i,:),size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages);
        %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(:,i,:),size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages);
        %         %         end
        %         handles.Current.MergedInd = 1;
        if get(handles.MergeChannel1Check,'Value')
            mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
        else
            mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
        end
        if get(handles.MergeChannel2Check,'Value')
            mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
        else
            mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
        end
        handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
        axes(handles.MergedData);
        image(handles.Current.mergedImage);
        xlabel('Frame Number');
        %         axis image
        set(handles.SliceCount,'Visible','on');
        set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,2))]);
        set(handles.LocationDrawCheck,'Visible','on');
        %         ax(3) = [];
        %         linkaxes(ax);
    end
% end
if get(handles.LocationDrawCheck,'Value')
    if get(handles.MergeViewPop,'Value') == 2
        tyline = [0 size(handles.Current.Image1,2)];
        txline = [handles.Current.Ind handles.Current.Ind];
        xline = [0 size(handles.Current.Image1,2)];
        yline = [handles.Current.MergedInd handles.Current.MergedInd];
    elseif get(handles.MergeViewPop,'Value') == 3
        tyline = [0 size(handles.Current.Image1,1)];
        txline = [handles.Current.Ind handles.Current.Ind];
        yline = [0 size(handles.Current.Image1,1)];
        xline = [handles.Current.MergedInd handles.Current.MergedInd];
    end
    axes(handles.Channel1data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    axes(handles.Channel2data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    
    axes(handles.MergedData);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(txline,tyline,'w','LineWidth',0.75);
    hold off
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function WhiteSlider_C2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WhiteSlider_C2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function BlackSlider_C2_Callback(hObject, eventdata, handles)
% hObject    handle to BlackSlider_C2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
MaxBrightnessStack = get(handles.BlackSlider_C2,'Max');

BlackLevel_C2 = round(double((get(hObject,'Value')))*100)/100;
%Make sure we don't try to set the brightness outside of the slider limits
if BlackLevel_C2 < 0
    BlackLevel_C2 = 0;
elseif BlackLevel_C2 > MaxBrightnessStack
    BlackLevel_C2 = MaxBrightnessStack-0.1;
end
set(handles.BlackEdit_C2,'String',num2str(BlackLevel_C2));

WhiteLevel_C2 = str2double(get(handles.WhiteEdit_C2,'String'));

if BlackLevel_C2 >= WhiteLevel_C2
    WhiteLevel_C2 = BlackLevel_C2+0.1;
    set(handles.WhiteEdit_C2,'String',num2str(round(double(WhiteLevel_C2)*100)/100));
end

set(handles.BlackSlider_C2,'Val',BlackLevel_C2);
set(handles.WhiteSlider_C2,'Val',WhiteLevel_C2);
handles.Current.clims2 = [BlackLevel_C2 WhiteLevel_C2];


set(handles.Channel2data,'NextPlot','replacechildren');
axes(handles.Channel2data);
imagesc(handles.Current.Image2, handles.Current.clims2);       % plot image
% axis image;
%Plot ROI

if get(handles.showParticlesCheck, 'Value')
    if handles.isFitPSF2
        plotParticle(handles.Tracking2.Particles, ...
            handles.Current.Ind, handles.isFitPSF2);
    else
        plotParticle(handles.Tracking2.Centroids, ...
            handles.Current.Ind, handles.isFitPSF2);
    end
end

if get(handles.showTracksCheck, 'Value')
    plotTracks(handles.Current.Tracks2, handles.Current.Ind);
end

if isfield(handles.Process2,'ROIpos') && get(handles.showROIsCheck,'Value')
    plotROI(handles.Process2.ROIpos);
    
end
%Show merged data
% if isfield(handles,'MergedInd')

    if get(handles.MergeViewPop,'Value') == 1
        if get(handles.MergeChannel1Check,'Value')
            mergeChannel1 = handles.Current.Image1;
        else
            mergeChannel1 = zeros(size(handles.Current.Image1));
        end
        if get(handles.MergeChannel2Check,'Value')
            mergeChannel2 = handles.Current.Image2;
        else
            mergeChannel2 = zeros(size(handles.Current.Image2));
        end
        handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
        axes(handles.MergedData);
        image(handles.Current.mergedImage);
        axis image

    elseif get(handles.MergeViewPop,'Value') == 2 %X*T
        %         set(handles.mergeSlider,'Visible','on');
        %         set(handles.mergeSlider,'Enable','on');
        %
        %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,1));
        %
        %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,1));
        %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,1)-1) 1/(size(handles.Data1.imageStack(1).data,1)-1)]);
        %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[2 3 1]);
        %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[2 3 1]);
        %
        %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,1));
        %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,1));
        %         %         for i = 1:size(handles.Data2.imageStack(1).data,1)
        %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(i,:,:),size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages);
        %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(i,:,:),size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages);
        %         %         end
        %         handles.Current.MergedInd = 1;
        if get(handles.MergeChannel1Check,'Value')
            mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
        else
            mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
        end
        if get(handles.MergeChannel2Check,'Value')
            mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
        else
            mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
        end
        handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
        axes(handles.MergedData);
        image(handles.Current.mergedImage);
        xlabel('Frame Number');
        %         axis image
        set(handles.SliceCount,'Visible','on');
        set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,1))]);
        set(handles.LocationDrawCheck,'Visible','on');
        %         ax(3) = [];
        %         linkaxes(ax);

    elseif get(handles.MergeViewPop,'Value') == 3 %Y*T
        %         set(handles.mergeSlider,'Visible','on');
        %         set(handles.mergeSlider,'Enable','on');
        %
        %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,2))
        %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,2));
        %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,2)-1) 1/(size(handles.Data1.imageStack(1).data,2)-1)]);
        %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[1 3 2]);
        %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[1 3 2]);
        %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,2));
        %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,2));
        %         %         for i = 1:size(handles.Data2.imageStack(1).data,2)
        %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(:,i,:),size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages);
        %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(:,i,:),size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages);
        %         %         end
        %         handles.Current.MergedInd = 1;
        if get(handles.MergeChannel1Check,'Value')
            mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
        else
            mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
        end
        if get(handles.MergeChannel2Check,'Value')
            mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
        else
            mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
        end
        handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
        axes(handles.MergedData);
        image(handles.Current.mergedImage);
        xlabel('Frame Number');
        %         axis image
        set(handles.SliceCount,'Visible','on');
        set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,2))]);
        set(handles.LocationDrawCheck,'Visible','on');
        %         ax(3) = [];
        %         linkaxes(ax);
    end
% end
if get(handles.LocationDrawCheck,'Value')
    if get(handles.MergeViewPop,'Value') == 2
        tyline = [0 size(handles.Current.Image1,2)];
        txline = [handles.Current.Ind handles.Current.Ind];
        xline = [0 size(handles.Current.Image1,2)];
        yline = [handles.Current.MergedInd handles.Current.MergedInd];
    elseif get(handles.MergeViewPop,'Value') == 3
        tyline = [0 size(handles.Current.Image1,1)];
        txline = [handles.Current.Ind handles.Current.Ind];
        yline = [0 size(handles.Current.Image1,1)];
        xline = [handles.Current.MergedInd handles.Current.MergedInd];
    end
    axes(handles.Channel1data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    axes(handles.Channel2data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    
    axes(handles.MergedData);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(txline,tyline,'w','LineWidth',0.75);
    hold off
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function BlackSlider_C2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BlackSlider_C2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function BlackEdit_C2_Callback(hObject, eventdata, handles)
% hObject    handle to BlackEdit_C2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BlackEdit_C2 as text
%        str2double(get(hObject,'String')) returns contents of BlackEdit_C2 as a double
MaxBrightnessStack = get(handles.BlackSlider_C2,'Max');

BlackLevel_C2 = str2double(get(hObject,'String'));
%Make sure we don't try to set the brightness outside of the slider limits
if BlackLevel_C2 < 0
    BlackLevel_C2 = 0;
elseif BlackLevel_C2 > MaxBrightnessStack
    BlackLevel_C2 = MaxBrightnessStack-0.1;
end

set(hObject,'String',num2str(round(double(BlackLevel_C2)*100)/100));
WhiteLevel_C2 = str2double(get(handles.WhiteEdit_C2,'String'));

if BlackLevel_C2 >= WhiteLevel_C2
    WhiteLevel_C2 = BlackLevel_C2+0.1;
    set(handles.WhiteEdit_C2,'String',num2str(round(double(WhiteLevel_C2)*100)/100));
end

set(handles.BlackSlider_C2,'Val',BlackLevel_C2);
set(handles.WhiteSlider_C2,'Val',WhiteLevel_C2);
handles.Current.clims2 = [BlackLevel_C2 WhiteLevel_C2];


set(handles.Channel2data,'NextPlot','replacechildren');
axes(handles.Channel2data);
imagesc(handles.Current.Image2, handles.Current.clims2);       % plot image
% axis image;
%Plot ROI

if get(handles.showParticlesCheck, 'Value')
    if handles.isFitPSF2
        plotParticle(handles.Tracking2.Particles, ...
            handles.Current.Ind, handles.isFitPSF2);
    else
        plotParticle(handles.Tracking2.Centroids, ...
            handles.Current.Ind, handles.isFitPSF2);
    end
end

if get(handles.showTracksCheck, 'Value')
    plotTracks(handles.Current.Tracks2, handles.Current.Ind);
end

if isfield(handles.Process2,'ROIpos') && get(handles.showROIsCheck,'Value')
    plotROI(handles.Process2.ROIpos);
    
end


if get(handles.MergeViewPop,'Value') == 1
    if get(handles.MergeChannel1Check,'Value')
        mergeChannel1 = handles.Current.Image1;
    else
        mergeChannel1 = zeros(size(handles.Current.Image1));
    end
    if get(handles.MergeChannel2Check,'Value')
        mergeChannel2 = handles.Current.Image2;
    else
        mergeChannel2 = zeros(size(handles.Current.Image2));
    end
    handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    axis image
elseif get(handles.MergeViewPop,'Value') == 2 %X*T
    %         set(handles.mergeSlider,'Visible','on');
    %         set(handles.mergeSlider,'Enable','on');
    %
    %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,1));
    %
    %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,1));
    %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,1)-1) 1/(size(handles.Data1.imageStack(1).data,1)-1)]);
    %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[2 3 1]);
    %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[2 3 1]);
    %
    %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,1));
    %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,1));
    %         %         for i = 1:size(handles.Data2.imageStack(1).data,1)
    %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(i,:,:),size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages);
    %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(i,:,:),size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages);
    %         %         end
    %         handles.Current.MergedInd = 1;
    if get(handles.MergeChannel1Check,'Value')
        mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
    else
        mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
    end
    if get(handles.MergeChannel2Check,'Value')
        mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
    else
        mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
    end
    handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    xlabel('Frame Number');
    %         axis image
    set(handles.SliceCount,'Visible','on');
    set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,1))]);
    set(handles.LocationDrawCheck,'Visible','on');
    %         ax(3) = [];
    %         linkaxes(ax);
    
elseif get(handles.MergeViewPop,'Value') == 3 %Y*T
    %         set(handles.mergeSlider,'Visible','on');
    %         set(handles.mergeSlider,'Enable','on');
    %
    %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,2))
    %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,2));
    %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,2)-1) 1/(size(handles.Data1.imageStack(1).data,2)-1)]);
    %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[1 3 2]);
    %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[1 3 2]);
    %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,2));
    %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,2));
    %         %         for i = 1:size(handles.Data2.imageStack(1).data,2)
    %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(:,i,:),size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages);
    %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(:,i,:),size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages);
    %         %         end
    %         handles.Current.MergedInd = 1;
    if get(handles.MergeChannel1Check,'Value')
        mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
    else
        mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
    end
    if get(handles.MergeChannel2Check,'Value')
        mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
    else
        mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
    end
    handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    xlabel('Frame Number');
    %         axis image
    set(handles.SliceCount,'Visible','on');
    set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,2))]);
    set(handles.LocationDrawCheck,'Visible','on');
    %         ax(3) = [];
    %         linkaxes(ax);
end

if get(handles.LocationDrawCheck,'Value')
    if get(handles.MergeViewPop,'Value') == 2
        tyline = [0 size(handles.Current.Image1,2)];
        txline = [handles.Current.Ind handles.Current.Ind];
        xline = [0 size(handles.Current.Image1,2)];
        yline = [handles.Current.MergedInd handles.Current.MergedInd];
    elseif get(handles.MergeViewPop,'Value') == 3
        tyline = [0 size(handles.Current.Image1,1)];
        txline = [handles.Current.Ind handles.Current.Ind];
        yline = [0 size(handles.Current.Image1,1)];
        xline = [handles.Current.MergedInd handles.Current.MergedInd];
    end
    axes(handles.Channel1data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    axes(handles.Channel2data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    
    axes(handles.MergedData);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(txline,tyline,'w','LineWidth',0.75);
    hold off
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function BlackEdit_C2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BlackEdit_C2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function WhiteEdit_C2_Callback(hObject, eventdata, handles)
% hObject    handle to WhiteEdit_C2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WhiteEdit_C2 as text
%        str2double(get(hObject,'String')) returns contents of WhiteEdit_C2 as a double
MaxBrightnessStack = get(handles.BlackSlider_C2,'Max');

WhiteLevel_C2 = str2double(get(hObject,'String'));
%Make sure we don't try to set the brightness outside of the slider limits
if WhiteLevel_C2 <= 0
    WhiteLevel_C2 = 0.1;
elseif WhiteLevel_C2 > MaxBrightnessStack
    WhiteLevel_C2 = MaxBrightnessStack;
end

set(hObject,'String',num2str(round(double(WhiteLevel_C2)*100)/100));
BlackLevel_C2 = str2double(get(handles.BlackEdit_C2,'String'));

if BlackLevel_C2 >= WhiteLevel_C2
    BlackLevel_C2 = WhiteLevel_C2-0.1;
    set(handles.BlackEdit_C2,'String',num2str(round(double(BlackLevel_C2)*100)/100));
end

set(handles.BlackSlider_C2,'Val',BlackLevel_C2);
set(handles.WhiteSlider_C2,'Val',WhiteLevel_C2);
handles.Current.clims2 = [BlackLevel_C2 WhiteLevel_C2];


set(handles.Channel2data,'NextPlot','replacechildren');
axes(handles.Channel2data);
imagesc(handles.Current.Image2, handles.Current.clims2);       % plot image
% axis image;
%Plot ROI

if get(handles.showParticlesCheck, 'Value')
    if handles.isFitPSF2
        plotParticle(handles.Tracking2.Particles, ...
            handles.Current.Ind, handles.isFitPSF2);
    else
        plotParticle(handles.Tracking2.Centroids, ...
            handles.Current.Ind, handles.isFitPSF2);
    end
end

if get(handles.showTracksCheck, 'Value')
    plotTracks(handles.Current.Tracks2, handles.Current.Ind);
end

if isfield(handles.Process2,'ROIpos') && get(handles.showROIsCheck,'Value')
    plotROI(handles.Process2.ROIpos);
    
end

if get(handles.MergeViewPop,'Value') == 1
    if get(handles.MergeChannel1Check,'Value')
        mergeChannel1 = handles.Current.Image1;
    else
        mergeChannel1 = zeros(size(handles.Current.Image1));
    end
    if get(handles.MergeChannel2Check,'Value')
        mergeChannel2 = handles.Current.Image2;
    else
        mergeChannel2 = zeros(size(handles.Current.Image2));
    end
    handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    axis image
elseif get(handles.MergeViewPop,'Value') == 2 %X*T
    %         set(handles.mergeSlider,'Visible','on');
    %         set(handles.mergeSlider,'Enable','on');
    %
    %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,1));
    %
    %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,1));
    %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,1)-1) 1/(size(handles.Data1.imageStack(1).data,1)-1)]);
    %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[2 3 1]);
    %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[2 3 1]);
    %
    %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,1));
    %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,1));
    %         %         for i = 1:size(handles.Data2.imageStack(1).data,1)
    %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(i,:,:),size(handles.Data1.imageStack(1).data,2),handles.Data1.nImages);
    %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(i,:,:),size(handles.Data2.imageStack(1).data,2),handles.Data2.nImages);
    %         %         end
    %         handles.Current.MergedInd = 1;
    if get(handles.MergeChannel1Check,'Value')
        mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
    else
        mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
    end
    if get(handles.MergeChannel2Check,'Value')
        mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
    else
        mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
    end
    handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    xlabel('Frame Number');
    %         axis image
    set(handles.SliceCount,'Visible','on');
    set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,1))]);
    set(handles.LocationDrawCheck,'Visible','on');
    %         ax(3) = [];
    %         linkaxes(ax);
    
elseif get(handles.MergeViewPop,'Value') == 3 %Y*T
    %         set(handles.mergeSlider,'Visible','on');
    %         set(handles.mergeSlider,'Enable','on');
    %
    %         set(handles.mergeSlider,'Min',1,'Max',size(handles.Data1.imageStack(1).data,2))
    %         set(handles.mergeSlider,'Value',size(handles.Data1.imageStack(1).data,2));
    %         set(handles.mergeSlider,'SliderStep',[1/(size(handles.Data1.imageStack(1).data,2)-1) 1/(size(handles.Data1.imageStack(1).data,2)-1)]);
    %         handles.Current.mergedStack.Image1 = permute(handles.mergedStack.Image1,[1 3 2]);
    %         handles.Current.mergedStack.Image2 = permute(handles.mergedStack.Image2,[1 3 2]);
    %         %         handles.Current.mergedStack.Image1 = zeros(size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages,size(handles.Data1.imageStack(1).data,2));
    %         %         handles.Current.mergedStack.Image2 = zeros(size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages,size(handles.Data2.imageStack(1).data,2));
    %         %         for i = 1:size(handles.Data2.imageStack(1).data,2)
    %         %             handles.Current.mergedStack.Image1(:,:,i) = reshape(handles.mergedStack.Image1(:,i,:),size(handles.Data1.imageStack(1).data,1),handles.Data1.nImages);
    %         %             handles.Current.mergedStack.Image2(:,:,i) = reshape(handles.mergedStack.Image2(:,i,:),size(handles.Data2.imageStack(1).data,1),handles.Data2.nImages);
    %         %         end
    %         handles.Current.MergedInd = 1;
    if get(handles.MergeChannel1Check,'Value')
        mergeChannel1 = handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd);
    else
        mergeChannel1 = zeros(size(handles.Current.mergedStack.Image1(:,:,handles.Current.MergedInd)));
    end
    if get(handles.MergeChannel2Check,'Value')
        mergeChannel2 = handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd);
    else
        mergeChannel2 = zeros(size(handles.Current.mergedStack.Image2(:,:,handles.Current.MergedInd)));
    end
    handles.Current.mergedImage = merge2colorData(mergeChannel1,mergeChannel2,[1 2],handles.Current.clims1,handles.Current.clims2);
    axes(handles.MergedData);
    image(handles.Current.mergedImage);
    xlabel('Frame Number');
    %         axis image
    set(handles.SliceCount,'Visible','on');
    set(handles.SliceCount,'String',['Slice ' num2str(handles.Current.MergedInd) '/' num2str(size(handles.mergedStack.Image1,2))]);
    set(handles.LocationDrawCheck,'Visible','on');
    %         ax(3) = [];
    %         linkaxes(ax);
end
if get(handles.LocationDrawCheck,'Value')
    if get(handles.MergeViewPop,'Value') == 2
        tyline = [0 size(handles.Current.Image1,2)];
        txline = [handles.Current.Ind handles.Current.Ind];
        xline = [0 size(handles.Current.Image1,2)];
        yline = [handles.Current.MergedInd handles.Current.MergedInd];
    elseif get(handles.MergeViewPop,'Value') == 3
        tyline = [0 size(handles.Current.Image1,1)];
        txline = [handles.Current.Ind handles.Current.Ind];
        yline = [0 size(handles.Current.Image1,1)];
        xline = [handles.Current.MergedInd handles.Current.MergedInd];
    end
    axes(handles.Channel1data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    axes(handles.Channel2data);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(xline,yline,'w','LineWidth',0.75);
    hold off
    
    axes(handles.MergedData);
    delete(findobj(gca,'LineWidth',0.75));
    hold on
    plot(txline,tyline,'w','LineWidth',0.75);
    hold off
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function WhiteEdit_C2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WhiteEdit_C2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function MergeCobinding_Callback(hObject, eventdata, handles)
% hObject    handle to MergeCobinding (see GCBO)
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
%         defPar(5) = str2double(get(handles.PxSizeText,'String'));
        
    else
        defPar(1) = handles.Parameters.Bound.ThreshL;
        defPar(2) = handles.Parameters.Bound.ThreshH;
        defPar(3) = handles.Parameters.Bound.minBoundFrames;
        if isfield(handles.Parameters,'Acquisition') && isfield(handles.Parameters.Acquisition,'frameTime')
            defPar(4) = handles.Parameters.Acquisition.frameTime;
        else
            defPar(4) = 0.2;
        end
        
    end
    defPar(5) = str2double(get(handles.PxSizeText,'String'));
    
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
        num2str(defPar(3)), '0', num2str(defPar(4)),num2str(defPar(5)),'0'};
    prompt = {'Maximum Jump between consecutive frames',...
        'Maximum end to end distance',...
        'Minimum length of bound tracks',...
        'Points in the histogram after log sampling (0 = No Log Sampling)',...
        'Frame Time', 'Pixel Size (um)'...
        'Maximum Frames to Analyze (0 = all)'};
    dlgtitle = 'Parameters for the Survival Time distribution';
    AnalysisParam = inputdlg(prompt,dlgtitle, 1, defaults);
    if isempty(AnalysisParam)
        return
    end
    
    dlgtitle2 = 'Parameters for 2 channel interactions';
    prompt2 = {'Maximum Separation distance',...
        'Minimum frames for interaction'};
    defaults2 = {num2str(defPar(1)), num2str(defPar(3))};
    
    AnalysisParam2 = inputdlg(prompt2,dlgtitle2,1,defaults2);
    if isempty(AnalysisParam2)
        return
    end
    
    AnalysisParam = [AnalysisParam; AnalysisParam2];
    %     end
    
    % Send parameters and file names to the residence time GUI
    TwoColorBinding_GUI(AnalysisFlag,PathName,FileNames, AnalysisParam);
end


% --------------------------------------------------------------------
function SaveMatNoImages_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMatNoImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
defName = [handles.Data1.pathName, handles.Data1.fileName];
SearchStr = '(.*)\.\w*';
defName = regexprep(defName, SearchStr, '$1');
FilterSpec = {'*.mat'};
[FileNameOut,PathNameOut] = uiputfile(FilterSpec,'Save Variables',defName);
if FileNameOut ~=0
    Version = 1.0;
    
    %     save([PathNameOut, FileNameOut], 'Version');
    Results.Data1.fileName = handles.Data1.fileName;
    Results.Data1.pathName = handles.Data1.pathName;
    Results.Data1.nImages = handles.Data1.nImages;
    Results.Data1.clims = handles.Data1.clims;
    Results.Data1.MaxBrightnessStack = handles.Data1.MaxBrightnessStack;
    Results.Data1.clims_orig = handles.Data1.clims_orig;
    Results.Data1.clims_filt = handles.Data1.clims_filt;
    Results.Data2.fileName = handles.Data2.fileName;
    Results.Data2.pathName = handles.Data2.pathName;
    Results.Data2.nImages = handles.Data2.nImages;
    Results.Data2.clims = handles.Data2.clims;
    Results.Data2.MaxBrightnessStack = handles.Data2.MaxBrightnessStack;
    Results.Data2.clims_orig = handles.Data2.clims_orig;
    Results.Data2.clims_filt = handles.Data2.clims_filt;
    %     Results.Current = handles.Current;
    Results.Parameters = handles.Parameters;
    Results.isAlign = handles.isAlign;
    Results.Process1.ROIClass = handles.Process1.ROIClass;
    Results.Process1.ROIimage = handles.Process1.ROIimage;
    Results.Process1.ROIpos = handles.Process1.ROIpos;
    Results.Process1.ROIlabel = handles.Process1.ROIlabel;
    Results.Process1.clims = handles.Process1.clims;
    Results.Process1.MaxBrightnessStack = handles.Process1.MaxBrightnessStack;
    Results.Process2.ROIClass = handles.Process2.ROIClass;
    Results.Process2.ROIimage = handles.Process2.ROIimage;
    Results.Process2.ROIpos = handles.Process2.ROIpos;
    Results.Process2.ROIlabel = handles.Process2.ROIlabel;
    Results.Process2.clims = handles.Process2.clims;
    Results.Process2.MaxBrightnessStack = handles.Process2.MaxBrightnessStack;
    Results.Tracking1 = handles.Tracking1;
    Results.Tracking2 = handles.Tracking2;
    %     Results.Analysis1 = handles.Analysis1;
    %     Results.Analysis2 = handles.Analysis2;
    Results.PreAnalysis1 = handles.PreAnalysis1;
    Results.PreAnalysis2 = handles.PreAnalysis2;
    Results.isFitPSF1 = handles.isFitPSF1;
    Results.isFitPSF2 = handles.isFitPSF2;
    %     Results.mergedStack = handles.mergedStack;
    save([PathNameOut, FileNameOut], 'Version','Results', '-v7.3');
    
end
