function varargout = preProcess_GUI(varargin)
% PREPROCESS_GUI M-file for preProcess_GUI.fig
%      PREPROCESS_GUI, by itself, creates a new PREPROCESS_GUI or raises the existing
%      singleton*.
%
%      H = PREPROCESS_GUI returns the handle to a new PREPROCESS_GUI or the handle to
%      the existing singleton*.
%
%      PREPROCESS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREPROCESS_GUI.M with the given input arguments.
%
%      PREPROCESS_GUI('Property','Value',...) creates a new PREPROCESS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before preProcess_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to preProcess_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help preProcess_GUI

% Last Modified by GUIDE v2.5 24-Jul-2015 12:49:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @preProcess_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @preProcess_GUI_OutputFcn, ...
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


% --- Executes just before preProcess_GUI is made visible.
function preProcess_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to preProcess_GUI (see VARARGIN)

% Choose default command line output for preProcess_GUI
handles.output = hObject;



% Import data from the workspace
handles.Tracks = varargin{1};        % Import tracks.
handles.Stack = varargin{2};                % Import stacks of images.
handles.Centroids = varargin{3};     % Import Centroids
handles.pixelSize = varargin{4};     % Import pixelsize
handles.NFrames = varargin{5};       % Import number of frames.
handles.FileName = varargin{6};      %Import FileName.
handles.ROIpos = varargin{7};
% Find Axes
% handles.TrackAxes = findobj ('Tag', 'Tracks');
% handles.NPartAxes = findobj ('Tag', 'NParticles');
handles.Hist1Axes = findobj ('Tag', 'Hist1');
handles.Hist2Axes = findobj ('Tag', 'Hist2');


% Plot All the Tracks.
axes(handles.TracksAxes);
SumImage = zeros(length(handles.Stack(1).data(:,1)),...
    length(handles.Stack(1).data(1,:)));


for imageIx =  1:handles.NFrames;
    SumImage = SumImage +  double(handles.Stack(imageIx).data);
end;
handles.SumImage = SumImage;
colormap(gray);
imagesc(SumImage);
axis image
hold on;
hold all;

if size(handles.Tracks,2) < 8
    for i = 1:1:max(handles.Tracks(:,4))
        ix = find(handles.Tracks(:,4)==i);
        plot(handles.Tracks(ix,1),handles.Tracks(ix,2),'LineWidth',1);
        
    end;
    
    box on;
    title('Tracks of identified particles');
    xlabel('X-position [\mum]');
    ylabel('Y-position [\mum]');
    
    hold off;
    if size(handles.Centroids,2) > 12
        ROIindx = 13;
    elseif size(handles.Centroids,2) > 6 && size(handles.Centroids,2) < 13
        ROIindx = 7;
    else
        ROIindx = 0;
    end
    
    if ROIindx > 0
        nROIs = length(handles.ROIpos);
    else
        nROIs = 1;
    end
    
    % Calculate Number of particles for each frame
    axes(handles.NParticlesAxes);
    handles.Nparticles(:,2:nROIs+1) = CalculateNparticles(handles.Centroids, handles.NFrames,nROIs);
    handles.Nparticles(:,1) = 1:length(handles.Nparticles(:,2));
    plot(sum(handles.Nparticles(:,2:end),2), '+k');
    xlabel('Frame');
    ylabel('Number of particles');
    title('Number of detected particles per frame');
    
    % Add Airy Intensity to Tracks
    handles.Tracks = AddAiryIntensity(handles.Tracks, ...
        handles.pixelSize, handles.Stack);
    
    % Calculate histograms of intensities
    IntensityHist = [];
    [IntensityHist(:,2),IntensityHist(:,1)] = hist(handles.Tracks(:,6), 50:50:10000);
    [IntensityHist(:,4),IntensityHist(:,3)] = hist(handles.Tracks(:,7), 50:50:10000);
    [IntensityHist(:,6),IntensityHist(:,5)] = hist(handles.Tracks(:,6)-handles.Tracks(:,7), 50:50:10000);
    
    
    % Plot histograms of Intensity and Background
    
    axes(handles.Hist1)
    bar(IntensityHist(:,3),IntensityHist(:,4),'b','EdgeColor', 'none','FaceAlpha',0.5)
    hold on
    bar(IntensityHist(:,1),IntensityHist(:,2),'r', 'EdgeColor', 'none','FaceAlpha',0.5)
    hold off
    maxPlot = IntensityHist(find(IntensityHist(:,2) ~= 0, 1,'last'),1);
    xlim([0 maxPlot]);
    
    
    title({'Intensity and Background histograms', 'for tracked particles'});
    xlabel('Intensity [AU]')
    ylabel('Counts')
    legend('Background', 'Particle Intensity');
    
    
    axes(handles.Hist2)
    bar(IntensityHist(:,5),IntensityHist(:,6),'b')
    maxPlot = IntensityHist(find(IntensityHist(:,6) ~= 0, 1,'last'),5);
    xlim([0 maxPlot]);
    title({'Intensity (background subtracted) histograms', 'for tracked particles'});
    xlabel('Intensity [AU]')
    ylabel('Counts')
    
    
    handles.IntensityHist = IntensityHist;
    
    % Convert the Tracks to microns
    handles.Tracks_um = [];
    handles.Tracks_um(:,1:2) = handles.Tracks(:,1:2)*handles.pixelSize;
    handles.Tracks_um(:,3:7) = handles.Tracks(:,3:7);
    handles.Intensity = handles.Tracks_um(:,6:7);
else
    handles.Intensity = AddAiryIntensityDim(handles.Tracks,handles.pixelSize,handles.Stack);
    handles.Tracks_um = [];
    handles.Tracks_um(:,1:2) = handles.Tracks(:,1:2);
    handles.Tracks_um(:,3) = handles.Tracks(:,3)*handles.pixelSize;
    handles.Tracks_um(:,4:5) = handles.Tracks(:,4:5);
    handles.Tracks_um(:,6:8) = (handles.Tracks(:,6:8)-1)*handles.pixelSize;
    for i = 1:size(handles.Tracks_um,1)
        if handles.Tracks_um(i,4) == 1
            xbox1 = handles.Tracks(i,6);
            xbox2 = handles.Tracks(i,6) + handles.Tracks(i,3);
            ybox1 = handles.Tracks(i,7);
            ybox2 = handles.Tracks(i,8);
        else
            ybox1 = handles.Tracks(i,6);
            ybox2 = handles.Tracks(i,6) + handles.Tracks(i,3);
            xbox1 = handles.Tracks(i,7);
            xbox2 = handles.Tracks(i,8);
        end
        boxBound = [xbox1, ybox1; xbox1,ybox2; xbox2,ybox2;xbox2,ybox1;xbox1,ybox1];
        plot(boxBound(:,1),boxBound(:,2));
        
        
    end
    hold off;
    if size(handles.Centroids,2) > 12
        ROIindx = 13;
    elseif size(handles.Centroids,2) > 6 && size(handles.Centroids,2) < 13
        ROIindx = 7;
    else
        ROIindx = 0;
    end
    
    if ROIindx > 0
        nROIs = max(handles.Centroids(:,ROIindx));
    else
        nROIs = 1;
    end
    
    % Calculate Number of particles for each frame
    axes(handles.NParticlesAxes);
    handles.Nparticles(:,2:nROIs+1) = CalculateNparticles(handles.Centroids, handles.NFrames);
    handles.Nparticles(:,1) = 1:length(handles.Nparticles(:,2));
    plot(sum(handles.Nparticles(:,2:end),2), '+k');
    xlabel('Frame');
    ylabel('Number of particles');
    title('Number of detected particles per frame');
    
    % Calculate histograms of intensities
    IntensityHist = [];
    [IntensityHist(:,2),IntensityHist(:,1)] = hist(handles.Intensity(:,2), 50:50:10000);
    [IntensityHist(:,4),IntensityHist(:,3)] = hist(handles.Intensity(:,3), 50:50:10000);
    [IntensityHist(:,6),IntensityHist(:,5)] = hist(handles.Intensity(:,2)-handles.Intensity(:,3), 50:50:10000);
    
    
    % Plot histograms of Intensity and Background
    
    axes(handles.Hist1)
    bar(IntensityHist(:,3),IntensityHist(:,4),'b','EdgeColor', 'none')
    hold on
    bar(IntensityHist(:,1),IntensityHist(:,2),'r', 'EdgeColor', 'none')
    hold off
    maxPlot = IntensityHist(find(IntensityHist(:,2) ~= 0, 1,'last'),1);
    xlim([0 maxPlot]);
    
    
    title({'Intensity and Background histograms', 'for tracked particles'});
    xlabel('Intensity [AU]')
    ylabel('Counts')
    legend('Background', 'Particle Intensity');
    
    
    axes(handles.Hist2)
    bar(IntensityHist(:,5),IntensityHist(:,6),'b')
    maxPlot = IntensityHist(find(IntensityHist(:,6) ~= 0, 1,'last'),5);
    xlim([0 maxPlot]);
    title({'Intensity (background subtracted) histograms', 'for tracked particles'});
    xlabel('Intensity [AU]')
    ylabel('Counts')
    
    
    handles.IntensityHist = IntensityHist;
end

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes preProcess_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = preProcess_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.Tracks_um;
varargout{2} = handles.Nparticles;
varargout{3} = handles.IntensityHist;


% --- Executes on button press in MakeMovie.
function MakeMovie_Callback(hObject, eventdata, handles)
% hObject    handle to MakeMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure;
nFrames = handles.NFrames;
clim = [min(handles.Stack(end).data(:)),...
    max(handles.Stack(end).data(:))];

green = gray;
green(:,1) = 0;
green(:,3) = 0;


%%%Added by David Garcia
DefaultName = handles.FileName;
SearchStr = '(.*)\.\w*';
DefaultName = regexprep(DefaultName, SearchStr, '$1_Movie');
FilterSpec = {'*.avi'};
[FileNameOut,PathNameOut] = uiputfile(FilterSpec,'Save Movie',DefaultName);
v = VideoWriter(strcat(PathNameOut,FileNameOut),'Uncompressed AVI');
open(v)
    
for i = 1:nFrames
    
    
    imshow(handles.Stack(i).data,clim,'InitialMagnification', 200);
    
    axis off
    ImIx = find (handles.Tracks(:,3) == i);
    
        if ~isempty(ImIx)
        pIx = handles.Tracks(ImIx,4);           %find the corrisponding particle index;
        for j =pIx'
            plotIx = find(handles.Tracks(:,4) == j & handles.Tracks(:,3) <= i);
            hold on;
            plot(handles.Tracks(plotIx,1),handles.Tracks(plotIx,2),'r','LineWidth',1);
            hold off;
        end
       
    end
     F(i) = getframe();
     writeVideo(v,F(i).cdata);
end
close(v)
    
    
% Save the movie
% DefaultName = handles.FileName;
% SearchStr = '(.*)\.\w*';
% DefaultName = regexprep(DefaultName, SearchStr, '$1_Movie');
% FilterSpec = {'*.avi'};
% [FileNameOut,PathNameOut] = uiputfile(FilterSpec,'Save Movie',DefaultName);
%  if FileNameOut ~= 0;
% %     movie2avi(F, [PathNameOut,FileNameOut], 'fps', 20);
%     v = VideoWriter(FileNameOut,'Uncompressed AVI');
%     open(v)
%     
%     writeVideo(v,F.cdata)
%     
%  end
    





% --- Executes on button press in CopyTracks.
function CopyTracks_Callback(hObject, eventdata, handles)
% hObject    handle to CopyTracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

num2clip(handles.Tracks_um);
msgbox({'The Array containing the tracks';...
    'has been copied to the clipboard'});


% --- Executes on button press in CopyNParticles.
function CopyNParticles_Callback(hObject, eventdata, handles)
% hObject    handle to CopyNParticles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


num2clip(handles.Nparticles);
msgbox({'The Array containing the Number of particles';...
    'per frame has been copied to the clipboard'});


% --- Executes on button press in CopyHistogram.
function CopyHistogram_Callback(hObject, eventdata, handles)
% hObject    handle to CopyHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

num2clip(handles.IntensityHist);
msgbox({'The Array containing the histogram of intensities';...
    'has been copied to the clipboard'});


% --- Executes on button press in SaveFigures.
function SaveFigures_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFigures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

IntensityHist = handles.IntensityHist;
Tracks = handles.Tracks;
Nparticles = handles.Nparticles;
SumImage = handles.SumImage;



% Plot Tracks
scrsz = get(0,'ScreenSize');
h1 = figure('Position',[.15*scrsz(3) .15*scrsz(4) .4*scrsz(3) .5*scrsz(4)]) ;

subplot(2,2,1)

colormap(gray);
imagesc(SumImage);
axis image
hold on;
hold all;


for i = 1:1:max(handles.Tracks(:,4))
    ix = find(handles.Tracks(:,4)==i);
    plot(handles.Tracks(ix,1),handles.Tracks(ix,2),'LineWidth',1);
 
end;

box on;
title('Tracks of identified particles');
xlabel('X-position [\mum]');
ylabel('Y-position [\mum]');

hold off;

subplot(2,2,2);

plot(sum(Nparticles(:,2:end),2), '+k');
xlabel('Frame');
ylabel('Number of particles');
title('Number of detected particles per frame');

subplot(2,2,3);
bar(IntensityHist(:,3),IntensityHist(:,4),'b','EdgeColor', 'none','FaceAlpha',0.5)
hold on
bar(IntensityHist(:,1),IntensityHist(:,2),'r', 'EdgeColor', 'none','FaceAlpha',0.5)
hold off
maxPlot = IntensityHist(find(IntensityHist(:,2) ~= 0, 1,'last'),1);
xlim([0 maxPlot]);


title({'Intensity and Background histograms', 'for tracked particles'});
xlabel('Intensity [AU]')
ylabel('Counts')
legend('Background', 'Particle Intensity');


subplot(2,2,4);
bar(IntensityHist(:,5),IntensityHist(:,6),'b')
maxPlot = IntensityHist(find(IntensityHist(:,6) ~= 0, 1,'last'),5);
xlim([0 maxPlot]);
title({'Intensity (background subtracted) histograms', 'for tracked particles'});
xlabel('Intensity [AU]')
ylabel('Counts')


% Save the figures
DefaultName = handles.FileName;
SearchStr = '(.*)\.\w*';
DefaultName = regexprep(DefaultName, SearchStr, '$1_preprocess');
FilterSpec = {'*.fig';'*.ai';'*.bpm';'*.eps';'*.jpeg';'*.pdf';'*.tif'};
[FileNameOut,PathNameOut] = uiputfile(FilterSpec,'Save Figure',DefaultName);
if FileNameOut ~= 0;
    saveas(h1, [PathNameOut,FileNameOut]);
end



% --- Executes on button press in Done.
function Done_Callback(hObject, eventdata, handles)
% hObject    handle to Done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close;
