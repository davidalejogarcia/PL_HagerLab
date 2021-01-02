function varargout = BoundMol_GUI(varargin)
% BOUNDMOL_GUI M-file for BoundMol_GUI.fig
%      BOUNDMOL_GUI, by itself, creates a new BOUNDMOL_GUI or raises the existing
%      singleton*.
%
%      H = BOUNDMOL_GUI returns the handle to a new BOUNDMOL_GUI or the handle to
%      the existing singleton*.
%
%      BOUNDMOL_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BOUNDMOL_GUI.M with the given input arguments.
%
%      BOUNDMOL_GUI('Property','Value',...) creates a new BOUNDMOL_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BoundMol_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BoundMol_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BoundMol_GUI

% Last Modified by GUIDE v2.5 24-Jan-2012 09:33:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BoundMol_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @BoundMol_GUI_OutputFcn, ...
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


% --- Executes just before BoundMol_GUI is made visible.
function BoundMol_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BoundMol_GUI (see VARARGIN)



% Import input
handles.Tracks = varargin{1};        % Import tracks.
handles.Parameters = varargin{2};    % Import Parameters.
handles.FrameTime = varargin{3};     % Import Frametime
handles.NParticles = varargin{4};    % Import NParticles
handles.NFrames = varargin{5};       % Import number of frames.
handles.FileName = varargin{6};      % Import FileName.
handles.ROIname = varargin{7};

% Identify Immobile Particles;
if size(handles.Tracks,2) < 8
    handles.Answer = questdlg({'Do you want to normalize the histogram';...
    'of  residence times to the bound fraction?'});

    [ImmTracks, Dummy] = calculateImmobileTracks...
        (handles.Tracks, handles.Parameters.ThreshL,...
        handles.Parameters.minBoundFrames, handles.Parameters.ThreshH, 0);


    % Plot immobile tracks;

    handles.BoundTracksP = findobj ('Tag', 'BoundTracks_plot');
    axes(handles.BoundTracksP);
    axis ij;
    hold on
    for i = 1:max(handles.Tracks(:,4));
        idx = find(handles.Tracks(:,4)==i);
        plot(handles.Tracks(idx,1), handles.Tracks(idx,2),'-ok', 'MarkerSize', 5);
    end



    if  ~isempty(ImmTracks);
        for i = 1:max(ImmTracks(:,4));
            idx = find(ImmTracks(:,4)==i);
            plot(ImmTracks(idx,1), ImmTracks(idx,2),'-og', 'MarkerSize', 5);
        end

    end

    title({'Identified Bound segments for the following settings';...
        ['Max Jump between consecutive frames = ',...
        num2str(handles.Parameters.ThreshL), '\mum'];
        ['Minimum number of consecutive frames = ',...
        num2str(handles.Parameters.minBoundFrames)];...
        ['Max Jump within ', num2str(handles.Parameters.minBoundFrames),...
        ' frames = ', num2str(handles.Parameters.ThreshH), '\mum']});

    xlabel('x [\mum]');
    ylabel('y [\mum]');
    box on;
    hold off
    figname = get(handles.figure1,'Name');
    figname = [figname ' ' handles.ROIname];
    set(handles.figure1,'Name',figname);


    % Calculate bound fraction;
    BoundFraction(:,1) = 1:handles.NFrames;
    if isempty(ImmTracks)
        BoundFraction(:,2) = 0;
        TrackLengthHist = [0 0];
    else
        BoundFraction(:,2) = zeros(1,handles.NFrames);
        % CountNumber of bound molecules per frame;
        for i = 1: handles.NFrames;
            BoundFraction(i, 2) = sum(ImmTracks(:,3) == i);
        end;


        BoundFraction(:,2) = BoundFraction(:,2)./handles.NParticles(:,2);
        idx = find(BoundFraction(:,2) == Inf);
        BoundFraction(idx, 2) = 1;



      % Calculate histogram of bound particles (cumulative)
        TrackLengthHist = ...
            calculateTrackLength(ImmTracks, handles.FrameTime,handles.Parameters.minBoundFrames);

      % Reduce the cumulative histogram to the exact histogram of residence
      % times
        TrackLengthHist(1:end-1,2) = ...
            (TrackLengthHist(1:end-1,2) - TrackLengthHist(2:end,2));

        switch handles.Answer
            case 'Yes'
        % Normalize the histogram of the residence times to the bound fraction;
        TrackLengthHist(:,2) = ...
            TrackLengthHist(:,2)*nanmean(BoundFraction(:,2))/sum(TrackLengthHist(:,2));
        end

    end

    handles.ImmTracks = ImmTracks;
    handles.BoundFraction = BoundFraction;
    handles.TrackLengthHist = TrackLengthHist;
    % Find Axes
    handles.BoundPlot = findobj ('Tag', 'BoundF_plot');
    handles.ResTPlot = findobj ('Tag', 'ResidenceTime_plot');

    %Plot Bound fraction;
    axes(handles.BoundPlot);
    plot(handles.BoundFraction(:,1), handles.BoundFraction(:,2));
    title({'Bound fraction over Time'; ...
        ['Average Bound Fraction = ', num2str(nanmean(handles.BoundFraction(:,2)))]});
    xlabel('Frame');
    ylabel('Bound Fraction');

    axes(handles.ResTPlot);
    bar(handles.TrackLengthHist(:,1), handles.TrackLengthHist(:,2));
    title('Histogram of residence times');
    xlabel('Residence Time (s)');
    switch handles.Answer
        case 'Yes'
            ylabel('Normalized Frequency');
        case 'No'
            ylabel('Counts');
    end
else
    ImmTracks = handles.Tracks(handles.Tracks(:,2) >= handles.Parameters.minBoundFrames,:);
    ImmTracks = ImmTracks(ImmTracks(:,3) <= handles.Parameters.ThreshL,:);
    if ~isempty(ImmTracks)
        TrackLength = ImmTracks(:,2);
        LongestTrack = max(TrackLength);
        ShortestTrack = handles.Parameters.minBoundFrames;
        
        TrackLengthHist = zeros(LongestTrack-ShortestTrack + 1,2);
        TrackLengthHist(:,1) = (ShortestTrack: LongestTrack)* handles.FrameTime;
%         TrackLengthHist(:,2) = hist(TrackLength,(ShortestTrack:LongestTrack));
        for  i = 1:length(TrackLength);
            TrackLengthHist(1:TrackLength(i)- ShortestTrack + 1,2) = ...
                TrackLengthHist(1:TrackLength(i) - ShortestTrack + 1,2)+1;
        end;
        % Reduce the cumulative histogram to the exact histogram of residence
      % times
        TrackLengthHist(1:end-1,2) = ...
            (TrackLengthHist(1:end-1,2) - TrackLengthHist(2:end,2));
    else
        TrackLengthHist = [];
    end
     handles.ImmTracks = ImmTracks;
    handles.BoundFraction = [];
    handles.TrackLengthHist = TrackLengthHist;
    handles.ResTPlot = findobj ('Tag', 'ResidenceTime_plot');
    axes(handles.ResTPlot);
    if ~isempty(handles.TrackLengthHist)
        bar(handles.TrackLengthHist(:,1), handles.TrackLengthHist(:,2));
        title('Histogram of residence times');
        xlabel('Residence Time (s)');
        ylabel('Counts');
    else
        errordlg('No data fits the criteria. Adjust the Min. Frames, or Maximum Jump and try again','No data for histogram');
        
    end
    
   
end


% Choose default command line output for BoundMol_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BoundMol_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BoundMol_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.BoundFraction;
varargout{2} = handles.TrackLengthHist;
varargout{3} = handles.ImmTracks; 



% --- Executes on button press in Copy_BoundF.
function Copy_BoundF_Callback(hObject, eventdata, handles)
% hObject    handle to Copy_BoundF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

num2clip(handles.BoundFraction);
msgbox({'The Array containing the bound fraction';...
    ' has been copied to the clipboard'});



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

num2clip(handles.TrackLengthHist);
msgbox({'The Array containing the histogram of residence times';...
    ' has been copied to the clipboard'});

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Create the figure and save it.

scrsz = get(0,'ScreenSize');
h1 = figure('Position',[.15*scrsz(3) .15*scrsz(4) .6*scrsz(3) .3*scrsz(4)]) ;


% Plot with the identified bound tracks;
subplot(1,3,1);
axis ij;
hold on;

for i = 1:max(handles.Tracks(:,4));
    idx = find(handles.Tracks(:,4)==i);
    plot(handles.Tracks(idx,1), handles.Tracks(idx,2),'-ok', 'MarkerSize', 3);
end



if  ~isempty(handles.ImmTracks);
    for i = 1:max(handles.ImmTracks(:,4));
        idx = find(handles.ImmTracks(:,4)==i);
        plot(handles.ImmTracks(idx,1), handles.ImmTracks(idx,2),'-og', 'MarkerSize', 3);
    end
    
end

title({'Identified Bound segments for the following settings';...
    ['Max Jump between consecutive frames = ',...
    num2str(handles.Parameters.ThreshL), '\mum'];
    ['Minimum number of consecutive frames = ',...
    num2str(handles.Parameters.minBoundFrames)];...
    ['Max Jump within ', num2str(handles.Parameters.minBoundFrames),...
    ' frames = ', num2str(handles.Parameters.ThreshH), '\mum']});

xlabel('x [\mum]');
ylabel('y [\mum]');
box on;
hold off

%Plot Bound fraction;
subplot(1, 3, 2);
plot(handles.BoundFraction(:,1), handles.BoundFraction(:,2));
title({'Bound fraction over Time'; ...
    ['Average Bound Fraction = ', num2str(nanmean(handles.BoundFraction(:,2)))]});
xlabel('Frame');
ylabel('Bound Fraction');

% Plot Histogram of residence times;
subplot(1, 3, 3);
bar(handles.TrackLengthHist(:,1), handles.TrackLengthHist(:,2));
title('Histogram of residence times');
xlabel('Residence Time');
switch handles.Answer
    case 'Yes'
        ylabel('Normalized Frequency');
    case 'No'
        ylabel('Counts');
end

DefaultName = handles.FileName;
SearchStr = '(.*)\.\w*';
DefaultName = regexprep(DefaultName, SearchStr, ['$1_' handles.ROIname '_BoundTracks']);
FilterSpec = {'*.fig';'*.ai';'*.bpm';'*.eps';'*.jpeg';'*.pdf';'*.tif'};
[FileNameOut,PathNameOut] = uiputfile(FilterSpec,'Save Figure',DefaultName);
 if FileNameOut ~= 0;
    saveas(h1, [PathNameOut,FileNameOut]);
    end



% --- Executes on button press in Close_button.
function Close_button_Callback(hObject, eventdata, handles)
% hObject    handle to Close_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close();
