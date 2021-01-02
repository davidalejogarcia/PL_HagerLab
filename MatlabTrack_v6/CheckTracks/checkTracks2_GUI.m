function varargout = checkTracks2_GUI(varargin)
% CHECKTRACKS2_GUI M-file for checkTracks2_GUI.fig
%      CHECKTRACKS2_GUI, by itself, creates a new CHECKTRACKS2_GUI or raises the existing
%      singleton*.
%
%      H = CHECKTRACKS2_GUI returns the handle to a new CHECKTRACKS2_GUI or the handle to
%      the existing singleton*.
%
%      CHECKTRACKS2_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHECKTRACKS2_GUI.M with the given input arguments.
%
%      CHECKTRACKS2_GUI('Property','Value',...) creates a new CHECKTRACKS2_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before checkTracks2_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to checkTracks2_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help checkTracks2_GUI

% Last Modified by GUIDE v2.5 15-Apr-2016 11:39:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @checkTracks2_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @checkTracks2_GUI_OutputFcn, ...
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




% --- Executes just before checkTracks2_GUI is made visible.
% ------------------------------------------------------------------------
function checkTracks2_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to checkTracks2_GUI (see VARARGIN)

% Choose default command line output for checkTracks2_GUI


handles.output = hObject;

% Import data from the workspace
handles.Stack = varargin{1};        % Import stack of images.
handles.Tracks = varargin{2};       % Import Tracks.
handles.Particles = varargin{3};    % Import Particles
handles.Process = varargin{4};
if ~isfield(handles.Process,'ROIpos')
    set(handles.roiShow,'Enable','off');
    set(handles.ROIselPop,'Enable','off');
else
    ROIstring{1,:} = 'All';
    for i = 1:length(handles.Process.ROIlabel)
        ROIstring{i+1,:} = handles.Process.ROIlabel{i,:};
    end
    set(handles.ROIselPop,'String',ROIstring);
end



handles.StackAxes = findobj ('Tag', 'Stack');   % Get handles to image axes
handles.TrackAxes = findobj ('Tag', 'Track');



% Get infos from the tracks

handles.nImages = size(handles.Stack,2);      % number of images
handles.nTracks = max(handles.Tracks(:,4));   % number of tracks
handles.outTracks = cell(1, handles.nTracks); % pre-allocate cells for
% output tracks
handles.outParticles = handles.Particles;      % copy particles to the output



for i = 1 : handles.nTracks
    idx = find(handles.Tracks(:,4) ==  i);
    handles.outTracks{i} = handles.Tracks(idx,:);
    
end

handles.currentTrIx = 1;
handles.currentTr = handles.outTracks{handles.currentTrIx};
% Set the track counter
set(handles.trackCounter, 'String',...
    ['Track ', num2str(handles.currentTrIx),'/',num2str(handles.nTracks)]);

% Set current Image
handles.currentIx = min(handles.currentTr(:,3));      % current image index
handles.currentImage = ...                            % current image
    handles.Stack(handles.currentIx).data;



for i = 1:handles.nImages
    
    handles.curImTimesWin(:,:,i) = handles.Stack(i).data;
end




pos1 = handles.currentTr(1,1:2);
t1 = handles.currentTr(1,3);
tend = handles.currentTr(end,3);
SpExtent = str2double(get(handles.SpaceWindow,'String'));
TiExtent = str2double(get(handles.TimeWindow,'String'));
xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));

twin = max(1,min(handles.currentIx,t1-TiExtent)):min(handles.nImages,max(handles.currentIx,tend+TiExtent));

ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
handles.ProjImage = ProjIm(:,twin);
% Set the image-counter
set(handles.imageCounter, 'String',...
    ['Image ', num2str(handles.currentIx),'/',num2str(handles.nImages)]);


% Prepare the image slider

set(handles.ImageSlider,'Min', 1);
set(handles.ImageSlider,'Max', handles.nImages);
set(handles.ImageSlider,'Value', handles.currentIx);
set(handles.ImageSlider, 'SliderStep', [1/(handles.nImages-1) 1/(handles.nImages-1)]);

% Deactivate the prev track button

set (handles.prevTrack, 'Enable', 'off');

%Determine Maximum intensity range of stack
MaxBrightness = zeros(handles.nImages,1);
for i = 1:handles.nImages
    MaxBrightness(i) = max(handles.Stack(i).data(:));
end
MaxBrightnessStack = 2*max(MaxBrightness); %give a 100% buffer above the maximum
set(handles.BlackValSlider,'Min',0,'Max',MaxBrightnessStack,'SliderStep',[1/(MaxBrightnessStack-1),10/(MaxBrightnessStack-1)]);
set(handles.WhiteValSlider,'Min',0,'Max',MaxBrightnessStack,'SliderStep',[1/(MaxBrightnessStack-1),10/(MaxBrightnessStack-1)]);


% Plot first frame of the track in the left plot

axes(handles.StackAxes);
colormap(gray);                                     % Set colormap
handles.clims = [min(min(handles.currentImage)) ... % Set color scale
    max(max(handles.currentImage))];

%handles.clims = [2000, 7000];
imagesc(handles.currentImage, handles.clims);       % plot image
axis image;

%update intensity slider positions & edit boxes
set(handles.BlackValSlider,'Value',handles.clims(1));
set(handles.WhiteValSlider,'Value',handles.clims(2));
set(handles.BlackValEdit,'String',num2str(round(double(handles.clims(1))*100)/100));
set(handles.WhiteValEdit,'String',num2str(round(double(handles.clims(2))*100)/100));

%Plot ROI
if get(handles.roiShow,'Value')
    plotROI(handles.Process.ROIpos,handles.StackAxes);
end
axes(handles.StackAxes);
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    hold on
    plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
    hold off
end


% Plot First frame + track overlay in the right plot

axes(handles.TrackAxes);
colormap(gray);                                     % Set colormap
% handles.clims = [min(min(handles.currentImage)) ... % Set color scale
% max(max(handles.currentImage))];

%handles.clims = [2000, 7000];
imagesc(twin,xwin,handles.ProjImage, handles.clims);
axis image;

plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
hold on;
plot([t1,t1],[xwin(1),xwin(end)],'g--');
plot([tend,tend],[xwin(1),xwin(end)],'r--');
plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
hold off;


% Update handles structure
guidata(hObject, handles);
outTracks = cell2mat(handles.outTracks');

setappdata(0,'outTracks',outTracks);
setappdata(0,'outParticles',handles.outParticles);

% UIWAIT makes checkTracks2_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% ------------------------------------------------------------------------



% --- Outputs from this function are returned to the command line.
function varargout = checkTracks2_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% % Get output from handles structure
varargout{1} = handles.output;
% varargout{1} = handles.out2{1};
% varargout{2} = handles.out2{2};
%
% % The figure can be deleted now
% close();





% --- Executes on slider movement.
% ------------------------------------------------------------------------
function ImageSlider_Callback(hObject, eventdata, handles)
% hObject    handle to ImageSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Set the current image

handles.currentIx = round(get(hObject,'Value'));
% Determine if we want to display Original or Filtered image
imageType = get(handles.imagePopup,'Value');

if imageType == 1 %Original
    handles.currentImage = ...                            % current image
        handles.Stack(handles.currentIx).data;
elseif imageType == 2 %Filtered
    handles.currentImage = ...
        handles.Process.filterStack(handles.currentIx).data;
end


% Set the image counter
set(handles.imageCounter, 'String',...
    ['Image ', num2str(handles.currentIx),'/',num2str(handles.nImages)]);



% Display current image

axes(handles.StackAxes);
imagesc(handles.currentImage, handles.clims);
axis image;

% If a particle is not present deactivate the 'Edit Track Point'
% button
if handles.currentIx > max(handles.currentTr(:,3))
    set(handles.StartTrackAtT,'Enable','off')
else
    set(handles.StartTrackAtT,'Enable','on')
end
if handles.currentIx < min(handles.currentTr(:,3))
    set(handles.EndTrackAtT,'Enable','off')
else
    set(handles.EndTrackAtT,'Enable','on')
end

idx = find (handles.currentTr(:,3) == handles.currentIx);
if isempty(idx)
    set(handles.editTrack, 'Enable', 'off');
else
    set(handles.editTrack, 'Enable', 'on');
end

%Plot the ROI
if get(handles.roiShow,'Value')
    plotROI(handles.Process.ROIpos,handles.StackAxes);
end


% Plot particle position



curProj = get(handles.ProjSelPop,'Value');

pos1 = handles.currentTr(1,1:2);
t1 = handles.currentTr(1,3);
tend = handles.currentTr(end,3);
SpExtent = str2double(get(handles.SpaceWindow,'String'));
TiExtent = str2double(get(handles.TimeWindow,'String'));
xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));


twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    hold on
    plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
    
    hold off
end
if curProj == 1
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,xwin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
    hold on;
    plot([t1,t1],[xwin(1),xwin(end)],'g--');
    plot([tend,tend],[xwin(1),xwin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
    
    
elseif curProj == 2
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
    hold on;
    plot([t1,t1],[ywin(1),ywin(end)],'g--');
    plot([tend,tend],[ywin(1),ywin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
elseif curProj == 3
    ProjIm = handles.currentImage;
    handles.ProjImage = ProjIm(ywin,xwin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(xwin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
    
end


% Update handles structure
guidata(hObject, handles);
%--------------------------------------------------------------------------


% --- Executes during object creation, after setting all properties.
function ImageSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in prevTrack.
% -------------------------------------------------------------------------

function prevTrack_Callback(hObject, eventdata, handles)
% hObject    handle to prevTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% Activate the Edit track point button
set(handles.editTrack, 'Enable', 'on');

% Set current track;

ROI_sel = get(handles.ROIselPop,'Value');
last_Tr = handles.currentTrIx;
handles.currentTrIx = handles.currentTrIx - 1;
handles.currentTr = handles.outTracks{handles.currentTrIx};
if ROI_sel > 1
    while(handles.currentTr(1,5)) ~= ROI_sel -1 && handles.currentTrIx > 1
        handles.currentTrIx = handles.currentTrIx - 1;
        handles.currentTr = handles.outTracks{handles.currentTrIx};
    end
    if handles.currentTrIx == 1 && handles.currentTr(1,5) ~= ROI_sel - 1
        msgbox('No more Tracks within the selected ROI', 'First Track');
        handles.currentTrIx = last_Tr;
        handles.currentTr = handles.outTracks{handles.currentTrIx};
        set(handles.prevTrack,'Enable','off');
    else
        handles.currentTrIxROI = handles.currentTrIxROI - 1;
        %         handles.nTracksROI = length(unique(Trk_tmp(:,4)));
        %
        %         set(handles.ROITrackCounter,'Visible','on');
        set(handles.ROITrackCounter,'String',['Track ', num2str(handles.currentTrIxROI), '/', num2str(handles.nTracksROI), ' in ROI']);
        
    end
end

% Activate the next. track button if moving from the last track
if handles.currentTrIx < handles.nTracks
    set(handles.nextTrack, 'Enable', 'On');
end

% Deactivate the prev button if moving to the first track
if handles.currentTrIx == 1
    set(handles.prevTrack, 'Enable', 'off');
end
% Set the track counter
set(handles.trackCounter, 'String',...
    ['Track ', num2str(handles.currentTrIx),'/',num2str(handles.nTracks)]);

% Set current Image
handles.currentIx = min(handles.currentTr(:,3));      % current image index

% Determine if we want to display Original or Filtered image
imageType = get(handles.imagePopup,'Value');

if imageType == 1 %Original
    handles.currentImage = ...                            % current image
        handles.Stack(handles.currentIx).data;
elseif imageType == 2 %Filtered
    handles.currentImage = ...
        handles.Process.filterStack(handles.currentIx).data;
end


% Set the image-counter
set(handles.imageCounter, 'String',...
    ['Image ', num2str(handles.currentIx),'/',num2str(handles.nImages)]);

% Set the image-slider
set(handles.ImageSlider,'Value', handles.currentIx);


% Plot first frame of the track in the left plot

axes(handles.StackAxes);
imagesc(handles.currentImage, handles.clims);       % plot image
axis image;

%Plot ROI
if get(handles.roiShow,'Value')
    plotROI(handles.Process.ROIpos,handles.StackAxes);
end




% Plot First frame + track overlay in the right plot

curProj = get(handles.ProjSelPop,'Value');

pos1 = handles.currentTr(1,1:2);
t1 = handles.currentTr(1,3);
tend = handles.currentTr(end,3);
SpExtent = str2double(get(handles.SpaceWindow,'String'));
TiExtent = str2double(get(handles.TimeWindow,'String'));
xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));


% twin = max(1,min(handles.currentIx,t1-TiExtent)):min(handles.nImages,max(handles.currentIx,tend+TiExtent));

twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    hold on
    plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
    
    hold off
end
if curProj == 1
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,xwin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
    hold on;
    plot([t1,t1],[xwin(1),xwin(end)],'g--');
    plot([tend,tend],[xwin(1),xwin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
    
    
elseif curProj == 2
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
    hold on;
    plot([t1,t1],[ywin(1),ywin(end)],'g--');
    plot([tend,tend],[ywin(1),ywin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
elseif curProj == 3
    ProjIm = handles.currentImage;
    handles.ProjImage = ProjIm(ywin,xwin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(xwin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
    
end

% Update handles structure
guidata(hObject, handles);

% -------------------------------------------------------------------------


% --- Executes on button press in nextTrack.
% ------------------------------------------------------------------------
function nextTrack_Callback(hObject, eventdata, handles)
% hObject    handle to nextTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% Activate the Edit track point button
set(handles.editTrack, 'Enable', 'on');

% Set curent track;
ROI_sel = get(handles.ROIselPop,'Value');
last_Tr = handles.currentTrIx;
handles.currentTrIx = handles.currentTrIx + 1;
handles.currentTr = handles.outTracks{handles.currentTrIx};
if ROI_sel > 1
    while(handles.currentTr(1,5)) ~= ROI_sel -1 && handles.currentTrIx < handles.nTracks
        handles.currentTrIx = handles.currentTrIx + 1;
        handles.currentTr = handles.outTracks{handles.currentTrIx};
    end
    if handles.currentTrIx == handles.nTracks && handles.currentTr(1,5) ~= ROI_sel - 1
        msgbox('No more Tracks within the selected ROI', 'Last Track');
        handles.currentTrIx = last_Tr;
        handles.currentTr = handles.outTracks{handles.currentTrIx};
        set(handles.nextTrack,'Enable','off');
    else
        handles.currentTrIxROI = handles.currentTrIxROI + 1;
        %         handles.nTracksROI = length(unique(Trk_tmp(:,4)));
        %
        %         set(handles.ROITrackCounter,'Visible','on');
        set(handles.ROITrackCounter,'String',['Track ', num2str(handles.currentTrIxROI), '/', num2str(handles.nTracksROI), ' in ROI']);
        
    end
end

% Activate the prev. track button if moving from track 1 to track 2
if handles.currentTrIx > 1
    set(handles.prevTrack, 'Enable', 'On');
end

% Deactivate the next button if moving to the last track
if handles.currentTrIx == handles.nTracks
    set(handles.nextTrack, 'Enable', 'off');
end

% Set the track counter
set(handles.trackCounter, 'String',...
    ['Track ', num2str(handles.currentTrIx),'/',num2str(handles.nTracks)]);

% Set current Image
handles.currentIx = min(handles.currentTr(:,3));      % current image index

% Determine if we want to display Original or Filtered image
imageType = get(handles.imagePopup,'Value');

if imageType == 1 %Original
    handles.currentImage = ...                            % current image
        handles.Stack(handles.currentIx).data;
elseif imageType == 2 %Filtered
    handles.currentImage = ...
        handles.Process.filterStack(handles.currentIx).data;
end


% Set the image-counter
set(handles.imageCounter, 'String',...
    ['Image ', num2str(handles.currentIx),'/',num2str(handles.nImages)]);

% Set the image-slider
set(handles.ImageSlider,'Value', handles.currentIx);


% Plot first frame of the track in the left plot

axes(handles.StackAxes);
imagesc(handles.currentImage, handles.clims);       % plot image
axis image;

%Plot ROI
if get(handles.roiShow,'Value')
    plotROI(handles.Process.ROIpos,handles.StackAxes);
end




% Plot First frame + track overlay in the right plot
curProj = get(handles.ProjSelPop,'Value');

pos1 = handles.currentTr(1,1:2);
t1 = handles.currentTr(1,3);
tend = handles.currentTr(end,3);
SpExtent = str2double(get(handles.SpaceWindow,'String'));
TiExtent = str2double(get(handles.TimeWindow,'String'));
xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));


twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    hold on
    plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
    
    hold off
end
if curProj == 1
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,xwin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
    hold on;
    plot([t1,t1],[xwin(1),xwin(end)],'g--');
    plot([tend,tend],[xwin(1),xwin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
    
    
elseif curProj == 2
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
    hold on;
    plot([t1,t1],[ywin(1),ywin(end)],'g--');
    plot([tend,tend],[ywin(1),ywin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
elseif curProj == 3
    ProjIm = handles.currentImage;
    handles.ProjImage = ProjIm(ywin,xwin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(xwin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
    
end


% Update handles structure
guidata(hObject, handles);

% -------------------------------------------------------------------------






% -------------------------------------------------------------------------
% --- Executes on button press in editTrack.


function editTrack_Callback(hObject, eventdata, handles)
% hObject    handle to editTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% check if Gaussian Fit is on

FitGauss = get(handles.gaussianFit, 'Value');


% edit track point position
set(handles.editTrack,'Enable', 'off')
set(handles.prevTrack,'Enable', 'off')
set(handles.nextTrack, 'Enable', 'off')
set(handles.deleteTrack, 'Enable', 'off')
set(handles.AdjCntr, 'Enable', 'off')
set(handles.done, 'Enable', 'off')
set(handles.addRemovePoint, 'Enable', 'off')
set(handles.gaussianFit, 'Enable', 'off')
set(handles.showParticle, 'Enable', 'off')
set(handles.JumpTo, 'Enable', 'off')
set(handles.ImageSlider,'Enable','off')

axes(handles.StackAxes)

hPoint = impoint(handles.StackAxes);
%     drawnow; pause(0.05);
newPos = getPosition(hPoint);
delete(hPoint);

% [newPos(1) newPos(2)] = ginput(1);
set(handles.editTrack,'Enable', 'on')
set(handles.prevTrack,'Enable', 'on')
set(handles.nextTrack, 'Enable', 'on')
set(handles.deleteTrack, 'Enable', 'on')
set(handles.AdjCntr, 'Enable', 'on')
set(handles.done, 'Enable', 'on')
set(handles.addRemovePoint, 'Enable', 'on')
set(handles.gaussianFit, 'Enable', 'on')
set(handles.showParticle, 'Enable', 'on')
set(handles.JumpTo, 'Enable', 'on')
set(handles.ImageSlider,'Enable','on')



if FitGauss % refine position with Gaussian fitting
    
    newPos = singlePeak_FIT(handles.Stack(handles.currentIx).data, newPos);
    
end

% update output track
idx = find(handles.currentTr(:,3) == handles.currentIx);
oldPos = handles.currentTr(idx, [1 2]);
handles.currentTr(idx,[1 2]) = newPos;
handles.outTracks{handles.currentTrIx} = handles.currentTr;

% update output particles
if length (handles.outParticles(1,:)) > 10
    PcoordX = 10;
    PcoordY = 11;
    ROIindx = 13;
else
    PcoordX = 1;
    PcoordY = 2;
    ROIindx = 7;
end
pIx = find(handles.outParticles(:,PcoordX) == oldPos(1) & ...
    handles.outParticles(:,PcoordY) == oldPos(2) & ...
    handles.outParticles (:,6) == handles.currentIx);

handles.outParticles(pIx,PcoordX) = newPos(1);
handles.outParticles(pIx,PcoordY) = newPos(2);
if isfield(handles.Process,'ROIimage')
    newPos = InsideROIcheck2([newPos handles.currentIx handles.currentTrIx],handles.Process.ROIimage);
    handles.outParticles(pIx,ROIindx) = newPos(5);
end

%Plot ROI
if get(handles.roiShow,'Value')
    plotROI(handles.Process.ROIpos,handles.StackAxes);
end

% plot new track
plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
curProj = get(handles.ProjSelPop,'Value');

pos1 = handles.currentTr(1,1:2);
t1 = handles.currentTr(1,3);
tend = handles.currentTr(end,3);
SpExtent = str2double(get(handles.SpaceWindow,'String'));
TiExtent = str2double(get(handles.TimeWindow,'String'));
xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));


twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    hold on
    plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
    
    hold off
end
if curProj == 1
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,xwin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
    hold on;
    plot([t1,t1],[xwin(1),xwin(end)],'g--');
    plot([tend,tend],[xwin(1),xwin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
    
    
elseif curProj == 2
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
    hold on;
    plot([t1,t1],[ywin(1),ywin(end)],'g--');
    plot([tend,tend],[ywin(1),ywin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
elseif curProj == 3
    ProjIm = handles.currentImage;
    handles.ProjImage = ProjIm(ywin,xwin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(xwin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
    
end

% Update handles structure
guidata(hObject, handles);
outTracks = cell2mat(handles.outTracks');

setappdata(0,'outTracks',outTracks);
setappdata(0,'outParticles',handles.outParticles);

% uiwait(handles.figure1);





%--------------------------------------------------------------------------
% --- Executes on button press in addRemovePoint.
function addRemovePoint_Callback(hObject, eventdata, handles)
% hObject    handle to addRemovePoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%find if there's a track element in the current frame
idx = find(handles.currentTr(:,3) == handles.currentIx);


if isempty(idx)
    % Add Particle if there's none
    set(handles.editTrack,'Enable', 'off')
    set(handles.prevTrack,'Enable', 'off')
    set(handles.nextTrack, 'Enable', 'off')
    set(handles.deleteTrack, 'Enable', 'off')
    set(handles.AdjCntr, 'Enable', 'off')
    set(handles.done, 'Enable', 'off')
    set(handles.addRemovePoint, 'Enable', 'off')
    set(handles.gaussianFit, 'Enable', 'off')
    set(handles.showParticle, 'Enable', 'off')
    set(handles.JumpTo, 'Enable', 'off')
    set(handles.ImageSlider,'Enable','off')
    
    hPoint = impoint(handles.StackAxes);
    newPos = getPosition(hPoint);
    %     drawnow; pause(0.05);
    delete(hPoint);
    
    %     [newPos(1) newPos(2)] = ginput(1);
    
    
    set(handles.editTrack,'Enable', 'on')
    set(handles.prevTrack,'Enable', 'on')
    set(handles.nextTrack, 'Enable', 'on')
    set(handles.deleteTrack, 'Enable', 'on')
    set(handles.AdjCntr, 'Enable', 'on')
    set(handles.done, 'Enable', 'on')
    set(handles.addRemovePoint, 'Enable', 'on')
    set(handles.gaussianFit, 'Enable', 'on')
    set(handles.showParticle, 'Enable', 'on')
    set(handles.JumpTo, 'Enable', 'on')
    set(handles.ImageSlider,'Enable','on')
    % check if Gaussian Fit is on
    
    FitGauss = get(handles.gaussianFit, 'Value');
    
    if FitGauss % refine position with Gaussian fitting
        
        [newPos,intensity,sigma,BG] = singlePeak_FIT(handles.Stack(handles.currentIx).data, newPos);
        
    end
    
    
    newPos(3) = handles.currentIx;
    newPos(4) = handles.currentTrIx;
    if isfield(handles.Process,'ROIimage')
        newPos = InsideROIcheck2(newPos,handles.Process.ROIimage);
    end
    NewParticle(1) = newPos(1);
    NewParticle(2) = newPos(2);
    if FitGauss
        NewParticle(3) = intensity;
        NewParticle(5) = intensity;
        
    else
        NewParticle(3) = handles.currentImage(round(newPos(2)),round(newPos(1)));
        NewParticle(5) = handles.currentImage(round(newPos(2)),round(newPos(1)));
    end
    NewParticle(6) = newPos(3);
    if size(handles.outParticles,2)> 10;
        if FitGauss
            NewParticle(7) = intensity;
            NewParticle(8) = sigma;
            NewParticle(9) = BG;
        end
        NewParticle(10) = newPos(1);
        NewParticle(11) = newPos(2);
        NewParticle(12) = 0;
        NewParticle(13) = newPos(5);
    else
        NewParticle(7) = newPos(5);
    end
    
    
    
    % check if the particle is added at the beginning of the track
    if handles.currentIx == min(handles.currentTr(:,3)) - 1;
        idxMerge = [];
        
        % check if there's another track with a particle in that
        % position [same frame, within one pixel]
        for iTrack = 1:handles.currentTrIx-1
            idxMerge = find(handles.outTracks{1,iTrack}(:,3) == handles.currentIx...
                & abs(handles.outTracks{1, iTrack}(:,1)-newPos(1))<1 ...
                & abs(handles.outTracks{1, iTrack}(:,2)-newPos(2))<1);
            
            if ~isempty(idxMerge)
                break
            end
        end
        
        % if there's already a particle at the selected position,
        %   ask if you want to merge the two tracks
        
        if ~isempty(idxMerge)
            
            trSegment1 = handles.outTracks{1, iTrack}(1:idxMerge,:);
            trSegment2 = handles.currentTr;
            trSegment2(:,4) = trSegment1(1,4);
            
            plotTwotrSegments(1,trSegment1, trSegment2, handles.TrackAxes,(4-get(handles.ProjSelPop,'Value')));
            button = questdlg(['Found another track (Track #',...
                num2str(trSegment1(1,4)),...
                ') at this position. Do you want to merge the two tracks?'],...
                'Merge Track','Yes','No','Yes');
            plotTwotrSegments(0);
            
            if strcmp(button, 'Yes')
                % Reindex particle identifier
                for i  = handles.currentTrIx+1:handles.nTracks
                    handles.outTracks{i}(:,4) =handles.outTracks{i}(:,4)-1;
                end
                
                % delete current track
                handles.outTracks(handles.currentTrIx) =  [];
                
                % update merged track
                handles.outTracks{iTrack} = [trSegment1; trSegment2];
                handles.currentTr = handles.outTracks{handles.currentTrIx};
                
                % update tracks number and track counter
                
                handles.nTracks = handles.nTracks-1;
                set(handles.trackCounter, 'String',...
                    ['Track ', num2str(handles.currentTrIx),'/',num2str(handles.nTracks)]);
            else
                PartExist = find(handles.outParticles(:,6) == handles.currentIx...
                    & abs(handles.outParticles(:,1)-newPos(1))<1 ...
                    & abs(handles.outParticles(:,2)-newPos(2))<1);
                handles.currentTr = [newPos;handles.currentTr];
                handles.outTracks{handles.currentTrIx} = handles.currentTr;
                if isempty(PartExist)
                    handles.outParticles = [handles.outParticles; NewParticle];
                end
                
            end
            
        else % if no merging is possible just add the particle
            PartExist = find(handles.outParticles(:,6) == handles.currentIx...
                & abs(handles.outParticles(:,1)-newPos(1))<1 ...
                & abs(handles.outParticles(:,2)-newPos(2))<1);
            handles.currentTr = [newPos;handles.currentTr];
            handles.outTracks{handles.currentTrIx} = handles.currentTr;
            if isempty(PartExist)
                handles.outParticles = [handles.outParticles; NewParticle];
            end
        end
        
        % check if the particle is added at the beginning of the track
        
    elseif handles.currentIx == max(handles.currentTr(:,3)) + 1;
        
        idxMerge = [];
        
        % check if there's another track with a particle in that
        % position [same frame, within one pixel]
        for iTrack = handles.currentTrIx+1: handles.nTracks
            idxMerge = find(handles.outTracks{1,iTrack}(:,3) == handles.currentIx...
                & abs(handles.outTracks{1, iTrack}(:,1)-newPos(1))<1 ...
                & abs(handles.outTracks{1, iTrack}(:,2)-newPos(2))<1);
            
            if ~isempty(idxMerge)
                break
            end
        end
        
        % if there's already a particle at the selected position,
        %   ask if you want to merge the two tracks
        
        if ~isempty(idxMerge)
            
            trSegment1 = handles.currentTr;
            trSegment2 = handles.outTracks{1, iTrack}(idxMerge:end,:);
            
            plotTwotrSegments(1,trSegment1, trSegment2, handles.TrackAxes,(4-get(handles.ProjSelPop,'Value')));
            button = questdlg(['Found another track (Track #',...
                num2str(trSegment2(1,4)),...
                ') at this position. Do you want to merge the two tracks?'],...
                'Merge Track','Yes','No','Yes');
            plotTwotrSegments(0);
            
            if strcmp(button, 'Yes')
                trSegment2(:,4) = trSegment1(1,4);
                % Reindex particle identifier
                for i  = iTrack + 1:handles.nTracks
                    handles.outTracks{i}(:,4) =handles.outTracks{i}(:,4)-1;
                end
                % delete track with index iTrack
                handles.outTracks(iTrack) = [];
                
                % update merged track
                handles.outTracks{handles.currentTrIx} = [trSegment1; trSegment2];
                handles.currentTr = handles.outTracks{handles.currentTrIx};
                
                % update tracks number and track counter
                
                handles.nTracks = handles.nTracks-1;
                set(handles.trackCounter, 'String',...
                    ['Track ', num2str(handles.currentTrIx),'/',num2str(handles.nTracks)]);
                
            end
            
            %  if there's no other particle at the selected position,
            %  just add the selected position
        else
            handles.currentTr = [handles.currentTr;newPos];
            handles.outTracks{handles.currentTrIx} = handles.currentTr;
            handles.outParticles = [handles.outParticles; NewParticle];
        end
        
    else
        msgbox('You can only add particles to frames close to the current track!');
    end
    
    %Plot ROI
    if get(handles.roiShow,'Value')
        plotROI(handles.Process.ROIpos,handles.StackAxes);
    end
    % Plot the current Track
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    curProj = get(handles.ProjSelPop,'Value');
    
    pos1 = handles.currentTr(1,1:2);
    t1 = handles.currentTr(1,3);
    tend = handles.currentTr(end,3);
    SpExtent = str2double(get(handles.SpaceWindow,'String'));
    TiExtent = str2double(get(handles.TimeWindow,'String'));
    xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
    ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));
    
    
    twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
    if get(handles.showParticle, 'Value')
        plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
        hold on
        plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
        
        hold off
    end
    if curProj == 1
        ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
        handles.ProjImage = ProjIm(:,twin);
        
        axes(handles.TrackAxes);
        colormap(gray);                                     % Set colormap
        
        imagesc(twin,xwin,handles.ProjImage, handles.clims);
        axis image;
        
        plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
        hold on;
        plot([t1,t1],[xwin(1),xwin(end)],'g--');
        plot([tend,tend],[xwin(1),xwin(end)],'r--');
        plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
        hold off;
        set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
        
        
    elseif curProj == 2
        ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
        handles.ProjImage = ProjIm(:,twin);
        axes(handles.TrackAxes);
        colormap(gray);                                     % Set colormap
        
        imagesc(twin,ywin,handles.ProjImage, handles.clims);
        axis image;
        
        plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
        hold on;
        plot([t1,t1],[ywin(1),ywin(end)],'g--');
        plot([tend,tend],[ywin(1),ywin(end)],'r--');
        plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
        hold off;
        set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
    elseif curProj == 3
        ProjIm = handles.currentImage;
        handles.ProjImage = ProjIm(ywin,xwin);
        axes(handles.TrackAxes);
        colormap(gray);                                     % Set colormap
        
        imagesc(xwin,ywin,handles.ProjImage, handles.clims);
        axis image;
        
        plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
        
    end
    
    
    % Update handles structure
    guidata(hObject, handles);
    
else
    % Remove particle if there's one
    
    % find if the particle to remove is at the edge of the track
    if handles.currentIx == min(handles.currentTr(:,3)) ...
            || handles.currentIx == max(handles.currentTr(:,3))
        
        
        % update output particles
        if length (handles.outParticles(1,:)) > 10
            PcoordX = 10;
            PcoordY = 11;
        else
            PcoordX = 1;
            PcoordY = 2;
        end
        oldPos = handles.currentTr(idx,[1 2]);
        pIx = find(handles.outParticles(:,PcoordX) == oldPos(1) & ...
            handles.outParticles(:,PcoordY) == oldPos(2) & ...
            handles.outParticles (:,6) == handles.currentIx);
        
        % if it is the case just remove the particle
        handles.currentTr(idx,:) = [];
        handles.outTracks{handles.currentTrIx} = handles.currentTr;
        
        
        
        handles.outParticles(pIx, :) = [];
        
        %Plot ROI
        if get(handles.roiShow,'Value')
            plotROI(handles.Process.ROIpos,handles.StackAxes);
        end
        plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
        curProj = get(handles.ProjSelPop,'Value');
        
        pos1 = handles.currentTr(1,1:2);
        t1 = handles.currentTr(1,3);
        tend = handles.currentTr(end,3);
        SpExtent = str2double(get(handles.SpaceWindow,'String'));
        TiExtent = str2double(get(handles.TimeWindow,'String'));
        xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
        ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));
        
        
        twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
        if get(handles.showParticle, 'Value')
            plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
            hold on
            plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
            
            hold off
        end
        if curProj == 1
            ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
            handles.ProjImage = ProjIm(:,twin);
            
            axes(handles.TrackAxes);
            colormap(gray);                                     % Set colormap
            
            imagesc(twin,xwin,handles.ProjImage, handles.clims);
            axis image;
            
            plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
            hold on;
            plot([t1,t1],[xwin(1),xwin(end)],'g--');
            plot([tend,tend],[xwin(1),xwin(end)],'r--');
            plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
            hold off;
            set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
            
            
        elseif curProj == 2
            ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
            handles.ProjImage = ProjIm(:,twin);
            axes(handles.TrackAxes);
            colormap(gray);                                     % Set colormap
            
            imagesc(twin,ywin,handles.ProjImage, handles.clims);
            axis image;
            
            plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
            hold on;
            plot([t1,t1],[ywin(1),ywin(end)],'g--');
            plot([tend,tend],[ywin(1),ywin(end)],'r--');
            plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
            hold off;
            set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
        elseif curProj == 3
            ProjIm = handles.currentImage;
            handles.ProjImage = ProjIm(ywin,xwin);
            axes(handles.TrackAxes);
            colormap(gray);                                     % Set colormap
            
            imagesc(xwin,ywin,handles.ProjImage, handles.clims);
            axis image;
            
            plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
            
        end
        
    else  % if the segment to cut is in the middle of the track Split the particle in two.
        
        trSegment1 = handles.currentTr(1:idx-1,:);
        trSegment2 = handles.currentTr(idx+1:end,:);
        trSegment2(:,4) = trSegment2(:,4) + 1;
        
        plotTwotrSegments(1,trSegment1, trSegment2, handles.TrackAxes,(4-get(handles.ProjSelPop,'Value')));
        button = questdlg('The track will be split in two - proceed?','Split Track','Yes','No','Yes');
        plotTwotrSegments(0);
        
        
        if strcmp(button, 'Yes')
            % update output particles
            if length (handles.outParticles(1,:)) > 10
                PcoordX = 10;
                PcoordY = 11;
            else
                PcoordX = 1;
                PcoordY = 2;
            end
            oldPos = handles.currentTr(idx,[1 2]);
            pIx = find(handles.outParticles(:,PcoordX) == oldPos(1) & ...
                handles.outParticles(:,PcoordY) == oldPos(2) & ...
                handles.outParticles (:,6) == handles.currentIx);
            handles.outParticles(pIx, :) = [];
            
            
            
            % Reindex particle identifier
            for i  = handles.currentTrIx+1:handles.nTracks
                handles.outTracks{i}(:,4) =handles.outTracks{i}(:,4)+1;
            end
            % concatenate tracks
            handles.outTracks = {handles.outTracks{1,1:handles.currentTrIx}, ...
                trSegment2, ...
                handles.outTracks{1,handles.currentTrIx+1:end}};
            
            % update current track
            handles.outTracks{handles.currentTrIx} = trSegment1;
            handles.currentTr = handles.outTracks{handles.currentTrIx};
            
            % update tracks number and track counter
            
            handles.nTracks = handles.nTracks+1;
            set(handles.trackCounter, 'String',...
                ['Track ', num2str(handles.currentTrIx),'/',num2str(handles.nTracks)]);
            
            
        end
        %Plot ROI
        if get(handles.roiShow,'Value')
            plotROI(handles.Process.ROIpos,handles.StackAxes);
        end
        % plot new track
        plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
        curProj = get(handles.ProjSelPop,'Value');
        
        pos1 = handles.currentTr(1,1:2);
        t1 = handles.currentTr(1,3);
        tend = handles.currentTr(end,3);
        SpExtent = str2double(get(handles.SpaceWindow,'String'));
        TiExtent = str2double(get(handles.TimeWindow,'String'));
        xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
        ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));
        
        
        twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
        if get(handles.showParticle, 'Value')
            plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
            hold on
            plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
            
            hold off
        end
        if curProj == 1
            ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
            handles.ProjImage = ProjIm(:,twin);
            
            axes(handles.TrackAxes);
            colormap(gray);                                     % Set colormap
            
            imagesc(twin,xwin,handles.ProjImage, handles.clims);
            axis image;
            
            plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
            hold on;
            plot([t1,t1],[xwin(1),xwin(end)],'g--');
            plot([tend,tend],[xwin(1),xwin(end)],'r--');
            plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
            hold off;
            set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
            
            
        elseif curProj == 2
            ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
            handles.ProjImage = ProjIm(:,twin);
            axes(handles.TrackAxes);
            colormap(gray);                                     % Set colormap
            
            imagesc(twin,ywin,handles.ProjImage, handles.clims);
            axis image;
            
            plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
            hold on;
            plot([t1,t1],[ywin(1),ywin(end)],'g--');
            plot([tend,tend],[ywin(1),ywin(end)],'r--');
            plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
            hold off;
            set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
        elseif curProj == 3
            ProjIm = handles.currentImage;
            handles.ProjImage = ProjIm(ywin,xwin);
            axes(handles.TrackAxes);
            colormap(gray);                                     % Set colormap
            
            imagesc(xwin,ywin,handles.ProjImage, handles.clims);
            axis image;
            
            plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
            
        end
        
        
        
        
        
        
        
        
    end
    % Update handles structure
    guidata(hObject, handles);
    outTracks = cell2mat(handles.outTracks');
    
    setappdata(0,'outTracks',outTracks);
    setappdata(0,'outParticles',handles.outParticles);
    
end
% uiwait(handles.figure1);







% --- Executes on button press in mergeTracks.
function mergeTracks_Callback(hObject, eventdata, handles)
% hObject    handle to mergeTracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in deleteTrack.
% -------------------------------------------------------------------------
function deleteTrack_Callback(hObject, eventdata, handles)
% hObject    handle to deleteTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


button = questdlg('Do you want to delete this track?','Delete Track','Yes','No','Yes');

if strcmp(button, 'Yes')
    
    % Reindex particle identifier
    for i  = handles.currentTrIx+1:handles.nTracks
        handles.outTracks{i}(:,4) =handles.outTracks{i}(:,4)-1;
    end
    
    % Delete Track
    handles.outTracks(handles.currentTrIx) = [];
    
    
    % Delete particles of the track
    % update output particles
    if length (handles.outParticles(1,:)) > 10
        PcoordX = 10;
        PcoordY = 11;
    else
        PcoordX = 1;
        PcoordY = 2;
    end
    
    % Define particles to be deleted
    indx = 1;
    for i = 1:length(handles.currentTr(:,1))
        i;
        if i == 46
            anig = 10;
        end
        if ~isempty(find(handles.outParticles(:,PcoordX) == handles.currentTr(i,1) & ...
                handles.outParticles(:,PcoordY) == handles.currentTr(i,2) & ...
                handles.outParticles (:,6) == handles.currentTr(i,3), 1))
            pIx(:,indx) = find(handles.outParticles(:,PcoordX) == handles.currentTr(i,1) & ...
                handles.outParticles(:,PcoordY) == handles.currentTr(i,2) & ...
                handles.outParticles (:,6) == handles.currentTr(i,3));
            indx = indx+1;
        end
        
    end
    pIx = pIx(:);
    pIx = pIx(pIx > 0);
    handles.outParticles(pIx, :) = [];
    
    
    
    %if you are at the last track update the index of the
    %current track
    if handles.currentTrIx == handles.nTracks
        handles.currentTrIx = handles.currentTrIx - 1;
    end
    handles.currentTr = handles.outTracks{handles.currentTrIx};
    handles.nTracks = handles.nTracks-1;
    
    % update  track counter
    
    set(handles.trackCounter, 'String',...
        ['Track ', num2str(handles.currentTrIx),'/',num2str(handles.nTracks)]);
    
    %Set the current frame to the first frame of the next track
    handles.currentIx = min(handles.currentTr(:,3));      % current image index
    % Determine if we want to display Original or Filtered image
    imageType = get(handles.imagePopup,'Value');
    
    if imageType == 1 %Original
        handles.currentImage = ...                            % current image
            handles.Stack(handles.currentIx).data;
    elseif imageType == 2 %Filtered
        handles.currentImage = ...
            handles.Process.filterStack(handles.currentIx).data;
    end
    
    
    % Set the image-counter
    set(handles.imageCounter, 'String',...
        ['Image ', num2str(handles.currentIx),'/',num2str(handles.nImages)]);
    
    % Set the image-slider
    set(handles.ImageSlider,'Value', handles.currentIx);
    
    % Plot first frame of the track in the left plot
    
    axes(handles.StackAxes);
    imagesc(handles.currentImage, handles.clims);       % plot image
    axis image;
    
    %Plot ROI
    if get(handles.roiShow,'Value')
        plotROI(handles.Process.ROIpos,handles.StackAxes);
    end
    % plot new track
    curProj = get(handles.ProjSelPop,'Value');
    
    pos1 = handles.currentTr(1,1:2);
    t1 = handles.currentTr(1,3);
    tend = handles.currentTr(end,3);
    SpExtent = str2double(get(handles.SpaceWindow,'String'));
    TiExtent = str2double(get(handles.TimeWindow,'String'));
    xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
    ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));
    
    
    twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
    if get(handles.showParticle, 'Value')
        plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
        hold on
        plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
        
        hold off
    end
    if curProj == 1
        ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
        handles.ProjImage = ProjIm(:,twin);
        
        axes(handles.TrackAxes);
        colormap(gray);                                     % Set colormap
        
        imagesc(twin,xwin,handles.ProjImage, handles.clims);
        axis image;
        
        plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
        hold on;
        plot([t1,t1],[xwin(1),xwin(end)],'g--');
        plot([tend,tend],[xwin(1),xwin(end)],'r--');
        plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
        hold off;
        set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
        
        
    elseif curProj == 2
        ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
        handles.ProjImage = ProjIm(:,twin);
        axes(handles.TrackAxes);
        colormap(gray);                                     % Set colormap
        
        imagesc(twin,ywin,handles.ProjImage, handles.clims);
        axis image;
        
        plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
        hold on;
        plot([t1,t1],[ywin(1),ywin(end)],'g--');
        plot([tend,tend],[ywin(1),ywin(end)],'r--');
        plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
        hold off;
        set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
    elseif curProj == 3
        ProjIm = handles.currentImage;
        handles.ProjImage = ProjIm(ywin,xwin);
        axes(handles.TrackAxes);
        colormap(gray);                                     % Set colormap
        
        imagesc(xwin,ywin,handles.ProjImage, handles.clims);
        axis image;
        
        plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
        
    end
    % Update handles structure
    guidata(hObject, handles);
    outTracks = cell2mat(handles.outTracks');
    
    setappdata(0,'outTracks',outTracks);
    setappdata(0,'outParticles',handles.outParticles);
    
end


% -------------------------------------------------------------------------




% --- Executes on button press in done.
% -------------------------------------------------------------------------
function done_Callback(hObject, eventdata, handles)

% merge the tracks together
% outTracks = cell2mat(handles.outTracks');

% outputtracks
% handles.out2{1} = outTracks;
% handles.out2{2} = handles.outParticles;

%update handles
% guidata(hObject, handles);

% setappdata(0,'outTracks',outTracks);
% setappdata(0,'outParticles',handles.outParticles);

close();
% Resume UI
% uiresume(handles.figure1);






% hObject    handle to done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in showParticle.
function showParticle_Callback(hObject, eventdata, handles)
% hObject    handle to showParticle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    pos1 = handles.currentTr(1,1:2);
    
    SpExtent = str2double(get(handles.SpaceWindow,'String'));
    TiExtent = str2double(get(handles.TimeWindow,'String'));
    xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
    ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));
    
    hold on
    plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
    hold off
else
    delete(findobj(handles.StackAxes,'Color','r'));
    delete(findobj(handles.StackAxes,'Color','w','LineStyle','--'));
end

% Hint: get(hObject,'Value') returns toggle state of showParticle


% --- Executes on button press in gaussianFit.
function gaussianFit_Callback(hObject, eventdata, handles)
% hObject    handle to gaussianFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gaussianFit


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)

% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

SliderMin = round(get(handles.ImageSlider,'Min'));
SliderMax = round(get(handles.ImageSlider,'Max'));
SliderIx = round(get(handles.ImageSlider,'Value'));

SliderStatus = get(handles.ImageSlider,'Enable');
%  Shortcut for NextFrame

if eventdata.Character == '.' & SliderIx < SliderMax  & strcmp(SliderStatus,'on')'
    set(handles.ImageSlider,'Value', SliderIx + 1);
    ImageSlider_Callback(handles.ImageSlider, eventdata, handles);
end

% Shortcut for PreviousFrame
if eventdata.Character == ','  & SliderIx > SliderMin & strcmp(SliderStatus,'on');
    set(handles.ImageSlider,'Value', SliderIx - 1);
    ImageSlider_Callback(handles.ImageSlider, eventdata, handles);
end

% Shortcut for Add/Remove Particle
if eventdata.Character == '/'
    addRemovePoint_Callback(handles.addRemovePoint, [], handles);
end

% Shortcut for edit Particle
if eventdata.Character =='m'
    
    idx = find (handles.currentTr(:,3) == handles.currentIx);
    if ~isempty(idx)
        editTrack_Callback(handles.editTrack, [], handles);
        
    end
end


% --- Executes on button press in JumpTo.
function JumpTo_Callback(hObject, eventdata, handles)
% hObject    handle to JumpTo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompt = {'Jump to track # :'};
dlg_title = 'Jump to track';
num_lines = 1;
def = {int2str(handles.currentIx)};
GoToTrack = inputdlg(prompt,dlg_title,num_lines,def);
GoToTrack =str2num(GoToTrack{1});

% Set curent track;
if GoToTrack > 0 & GoToTrack <= handles.nTracks
    handles.currentTrIx = GoToTrack;
    handles.currentTr = handles.outTracks{handles.currentTrIx};
    
    % activate/deactivate the Prev Track button
    if GoToTrack > 1
        set(handles.prevTrack, 'Enable', 'on');
        
        
    else
        set(handles.prevTrack, 'Enable', 'off');
    end
    
    % activate/deactivate the Next track button
    
    if GoToTrack < handles.nTracks
        set(handles.nextTrack, 'Enable', 'on');
        
        
    else
        set(handles.nextTrack, 'Enable', 'off');
    end
    
    % Set the track counter
    set(handles.trackCounter, 'String',...
        ['Track ', num2str(handles.currentTrIx),'/',num2str(handles.nTracks)]);
    
    % Set current Image
    handles.currentIx = min(handles.currentTr(:,3));      % current image index
    % Determine if we want to display Original or Filtered image
    imageType = get(handles.imagePopup,'Value');
    
    if imageType == 1 %Original
        handles.currentImage = ...                            % current image
            handles.Stack(handles.currentIx).data;
    elseif imageType == 2 %Filtered
        handles.currentImage = ...
            handles.Process.filterStack(handles.currentIx).data;
    end
    
    
    % Set the image-counter
    set(handles.imageCounter, 'String',...
        ['Image ', num2str(handles.currentIx),'/',num2str(handles.nImages)]);
    
    % Set the image-slider
    set(handles.ImageSlider,'Value', handles.currentIx);
    
    
    % Plot first frame of the track in the left plot
    
    axes(handles.StackAxes);
    imagesc(handles.currentImage, handles.clims);       % plot image
    axis image;
    %Plot ROI
    if get(handles.roiShow,'Value')
        plotROI(handles.Process.ROIpos,handles.StackAxes);
    end
    
    
    
    % Plot First frame + track overlay in the right plot
    
    curProj = get(handles.ProjSelPop,'Value');
    
    pos1 = handles.currentTr(1,1:2);
    t1 = handles.currentTr(1,3);
    tend = handles.currentTr(end,3);
    SpExtent = str2double(get(handles.SpaceWindow,'String'));
    TiExtent = str2double(get(handles.TimeWindow,'String'));
    xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
    ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));
    
    
    twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
    if get(handles.showParticle, 'Value')
        plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
        hold on
        plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
        
        hold off
    end
    if curProj == 1
        ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
        handles.ProjImage = ProjIm(:,twin);
        
        axes(handles.TrackAxes);
        colormap(gray);                                     % Set colormap
        
        imagesc(twin,xwin,handles.ProjImage, handles.clims);
        axis image;
        
        plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
        hold on;
        plot([t1,t1],[xwin(1),xwin(end)],'g--');
        plot([tend,tend],[xwin(1),xwin(end)],'r--');
        plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
        hold off;
        set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
        
        
    elseif curProj == 2
        ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
        handles.ProjImage = ProjIm(:,twin);
        axes(handles.TrackAxes);
        colormap(gray);                                     % Set colormap
        
        imagesc(twin,ywin,handles.ProjImage, handles.clims);
        axis image;
        
        plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
        hold on;
        plot([t1,t1],[ywin(1),ywin(end)],'g--');
        plot([tend,tend],[ywin(1),ywin(end)],'r--');
        plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
        hold off;
        set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
    elseif curProj == 3
        ProjIm = handles.currentImage;
        handles.ProjImage = ProjIm(ywin,xwin);
        axes(handles.TrackAxes);
        colormap(gray);                                     % Set colormap
        
        imagesc(xwin,ywin,handles.ProjImage, handles.clims);
        axis image;
        
        plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
        
    end
    
    
    % Update handles structure
    guidata(hObject, handles);
else
    h = errordlg(['Please select a proper track number']);
end;



% --- Executes on button press in AdjCntr.
function AdjCntr_Callback(hObject, eventdata, handles)

% Plot first frame of the track in the left plot

axes(handles.StackAxes);
colormap(gray);                                     % Set colormap
handles.clims = [min(min(handles.currentImage)) ... % Set color scale
    max(max(handles.currentImage))];

%update Intensity sliders & edit boxes
set(handles.BlackValSlider,'Val',handles.clims(1));
set(handles.WhiteValSlider,'Val',handles.clims(2));

set(handles.BlackValEdit,'String',num2str(round(double(handles.clims(1))*100)/100));
set(handles.WhiteValEdit,'String',num2str(round(double(handles.clims(2))*100)/100));


imagesc(handles.currentImage, handles.clims);       % plot image
axis image;
%Plot ROI
if get(handles.roiShow,'Value')
    plotROI(handles.Process.ROIpos,handles.StackAxes);
end




% Plot First frame + track overlay in the right plot

curProj = get(handles.ProjSelPop,'Value');

pos1 = handles.currentTr(1,1:2);
t1 = handles.currentTr(1,3);
tend = handles.currentTr(end,3);
SpExtent = str2double(get(handles.SpaceWindow,'String'));
TiExtent = str2double(get(handles.TimeWindow,'String'));
xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));



if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    hold on
    plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
    hold off
end
twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    hold on
    plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
    
    hold off
end
if curProj == 1
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,xwin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
    hold on;
    plot([t1,t1],[xwin(1),xwin(end)],'g--');
    plot([tend,tend],[xwin(1),xwin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
    
    
elseif curProj == 2
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
    hold on;
    plot([t1,t1],[ywin(1),ywin(end)],'g--');
    plot([tend,tend],[ywin(1),ywin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
elseif curProj == 3
    ProjIm = handles.currentImage;
    handles.ProjImage = ProjIm(ywin,xwin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(xwin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
    
end

% Update handles structure
guidata(hObject, handles);



% hObject    handle to AdjCntr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure



% outputtracks
% handles.out2{1} = outTracks;
% handles.out2{2} = handles.outParticles;


%update handles
% guidata(hObject, handles);

% merge the tracks together
% outTracks = cell2mat(handles.outTracks');
%
% setappdata(0,'outTracks',outTracks);
% setappdata(0,'outParticles',handles.outParticles);
%
% close();
% Resume UI
% uiresume(handles.figure1);


% --- Executes on selection change in imagePopup.
function imagePopup_Callback(hObject, eventdata, handles)
% hObject    handle to imagePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns imagePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imagePopup
% Determine if we want to display Original or Filtered image
imageType = get(hObject,'Value');
MaxBrightness = zeros(handles.nImages,1);

if imageType == 1 %Original
    handles.currentImage = ...                            % current image
        handles.Stack(handles.currentIx).data;
    for i = 1:handles.nImages
        MaxBrightness(i) = max(handles.Stack(i).data(:));
    end
    for i = 1:handles.nImages
        
        handles.curImTimesWin(:,:,i) = handles.Stack(i).data;
    end
    
elseif imageType == 2 %Filtered
    handles.currentImage = ...
        handles.Process.filterStack(handles.currentIx).data;
    for i = 1:handles.nImages
        MaxBrightness(i) = max(handles.Process.filterStack(i).data(:));
    end
    for i = 1:handles.nImages
        
        handles.curImTimesWin(:,:,i) = handles.Process.filterStack(i).data;
    end
    
end
MaxBrightnessStack = 2*max(MaxBrightness);

set(handles.BlackValSlider,'Min',0,'Max',MaxBrightnessStack,'SliderStep',[1/(MaxBrightnessStack-1),10/(MaxBrightnessStack-1)]);
set(handles.WhiteValSlider,'Min',0,'Max',MaxBrightnessStack,'SliderStep',[1/(MaxBrightnessStack-1),10/(MaxBrightnessStack-1)]);

axes(handles.StackAxes);
colormap(gray);                                     % Set colormap
handles.clims = [min(min(handles.currentImage)) ... % Set color scale
    max(max(handles.currentImage))];

%update intensity sliders & edit boxes
set(handles.BlackValSlider,'Value',handles.clims(1));
set(handles.WhiteValSlider,'Value',handles.clims(2));

set(handles.BlackValEdit,'String',num2str(handles.clims(1)));
set(handles.WhiteValEdit,'String',num2str(round(handles.clims(2)*100)/100));

imagesc(handles.currentImage, handles.clims);       % plot image
axis image;

%Plot ROI
if get(handles.roiShow,'Value')
    plotROI(handles.Process.ROIpos,handles.StackAxes);
end

% Plot particle position



curProj = get(handles.ProjSelPop,'Value');

pos1 = handles.currentTr(1,1:2);
t1 = handles.currentTr(1,3);
tend = handles.currentTr(end,3);
SpExtent = str2double(get(handles.SpaceWindow,'String'));
TiExtent = str2double(get(handles.TimeWindow,'String'));
xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));


% twin = max(1,min(handles.currentIx,t1-TiExtent)):min(handles.nImages,max(handles.currentIx,tend+TiExtent));
% if get(handles.showParticle, 'Value')
%     plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
%     hold on
%     plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
%     hold off
% end
twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    hold on
    plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
    
    hold off
end
if curProj == 1
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,xwin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
    hold on;
    plot([t1,t1],[xwin(1),xwin(end)],'g--');
    plot([tend,tend],[xwin(1),xwin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
    
    
elseif curProj == 2
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
    hold on;
    plot([t1,t1],[ywin(1),ywin(end)],'g--');
    plot([tend,tend],[ywin(1),ywin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
elseif curProj == 3
    ProjIm = handles.currentImage;
    handles.ProjImage = ProjIm(ywin,xwin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(xwin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
    
end

% Update handles structure
guidata(hObject, handles);


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
    set(handles.BlackValSlider,'Val',WhiteLevel);
end

set(handles.BlackValEdit,'String',num2str(round(double(BlackLevel)*100)/100));
set(handles.WhiteValEdit,'String',num2str(round(double(WhiteLevel)*100)/100));
handles.clims = [BlackLevel WhiteLevel];


axes(handles.StackAxes);
colormap(gray);                                     % Set colormap

imagesc(handles.currentImage, handles.clims);       % plot image
axis image;
%Plot ROI
if get(handles.roiShow,'Value')
    plotROI(handles.Process.ROIpos,handles.StackAxes);
end
% Plot particle position



curProj = get(handles.ProjSelPop,'Value');

pos1 = handles.currentTr(1,1:2);
t1 = handles.currentTr(1,3);
tend = handles.currentTr(end,3);
SpExtent = str2double(get(handles.SpaceWindow,'String'));
TiExtent = str2double(get(handles.TimeWindow,'String'));
xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));


% twin = max(1,min(handles.currentIx,t1-TiExtent)):min(handles.nImages,max(handles.currentIx,tend+TiExtent));
% if get(handles.showParticle, 'Value')
%     plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
%     hold on
%     plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
%     hold off
% end
twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    hold on
    plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
    
    hold off
end
if curProj == 1
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,xwin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
    hold on;
    plot([t1,t1],[xwin(1),xwin(end)],'g--');
    plot([tend,tend],[xwin(1),xwin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
    
    
elseif curProj == 2
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
    hold on;
    plot([t1,t1],[ywin(1),ywin(end)],'g--');
    plot([tend,tend],[ywin(1),ywin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
elseif curProj == 3
    ProjIm = handles.currentImage;
    handles.ProjImage = ProjIm(ywin,xwin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(xwin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
    
end
%update handles structure
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
handles.clims = [BlackLevel WhiteLevel];


axes(handles.StackAxes);
colormap(gray);                                     % Set colormap

imagesc(handles.currentImage, handles.clims);       % plot image
axis image;
%Plot ROI
if get(handles.roiShow,'Value')
    plotROI(handles.Process.ROIpos,handles.StackAxes);
end
% Plot particle position



curProj = get(handles.ProjSelPop,'Value');

pos1 = handles.currentTr(1,1:2);
t1 = handles.currentTr(1,3);
tend = handles.currentTr(end,3);
SpExtent = str2double(get(handles.SpaceWindow,'String'));
TiExtent = str2double(get(handles.TimeWindow,'String'));
xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));


% twin = max(1,min(handles.currentIx,t1-TiExtent)):min(handles.nImages,max(handles.currentIx,tend+TiExtent));
% if get(handles.showParticle, 'Value')
%     plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
%     hold on
%     plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
%     hold off
% end
twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    hold on
    plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
    
    hold off
end
if curProj == 1
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,xwin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
    hold on;
    plot([t1,t1],[xwin(1),xwin(end)],'g--');
    plot([tend,tend],[xwin(1),xwin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
    
    
elseif curProj == 2
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
    hold on;
    plot([t1,t1],[ywin(1),ywin(end)],'g--');
    plot([tend,tend],[ywin(1),ywin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
elseif curProj == 3
    ProjIm = handles.currentImage;
    handles.ProjImage = ProjIm(ywin,xwin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(xwin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
    
end
%update handles structure
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
handles.clims = [BlackLevel WhiteLevel];


axes(handles.StackAxes);
colormap(gray);                                     % Set colormap

imagesc(handles.currentImage, handles.clims);       % plot image
axis image;
%Plot ROI
if get(handles.roiShow,'Value')
    plotROI(handles.Process.ROIpos,handles.StackAxes);
end
% Plot particle position



curProj = get(handles.ProjSelPop,'Value');

pos1 = handles.currentTr(1,1:2);
t1 = handles.currentTr(1,3);
tend = handles.currentTr(end,3);
SpExtent = str2double(get(handles.SpaceWindow,'String'));
TiExtent = str2double(get(handles.TimeWindow,'String'));
xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));


% twin = max(1,min(handles.currentIx,t1-TiExtent)):min(handles.nImages,max(handles.currentIx,tend+TiExtent));
% if get(handles.showParticle, 'Value')
%     plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
%     hold on
%     plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
%     hold off
% end
twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    hold on
    plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
    
    hold off
end
if curProj == 1
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,xwin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
    hold on;
    plot([t1,t1],[xwin(1),xwin(end)],'g--');
    plot([tend,tend],[xwin(1),xwin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
    
    
elseif curProj == 2
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
    hold on;
    plot([t1,t1],[ywin(1),ywin(end)],'g--');
    plot([tend,tend],[ywin(1),ywin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
elseif curProj == 3
    ProjIm = handles.currentImage;
    handles.ProjImage = ProjIm(ywin,xwin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(xwin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
    
end

%update handles structure
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
    WhitekLevel = 0.1;
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
handles.clims = [BlackLevel WhiteLevel];


axes(handles.StackAxes);
colormap(gray);                                     % Set colormap

imagesc(handles.currentImage, handles.clims);       % plot image
axis image;
%Plot ROI
if get(handles.roiShow,'Value')
    plotROI(handles.Process.ROIpos,handles.StackAxes);
end
% Plot particle position



curProj = get(handles.ProjSelPop,'Value');

pos1 = handles.currentTr(1,1:2);
t1 = handles.currentTr(1,3);
tend = handles.currentTr(end,3);
SpExtent = str2double(get(handles.SpaceWindow,'String'));
TiExtent = str2double(get(handles.TimeWindow,'String'));
xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));


% twin = max(1,min(handles.currentIx,t1-TiExtent)):min(handles.nImages,max(handles.currentIx,tend+TiExtent));
% if get(handles.showParticle, 'Value')
%     plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
%     hold on
%     plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
%     hold off
% end
twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    hold on
    plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
    
    hold off
end
if curProj == 1
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,xwin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
    hold on;
    plot([t1,t1],[xwin(1),xwin(end)],'g--');
    plot([tend,tend],[xwin(1),xwin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
    
    
elseif curProj == 2
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
    hold on;
    plot([t1,t1],[ywin(1),ywin(end)],'g--');
    plot([tend,tend],[ywin(1),ywin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
elseif curProj == 3
    ProjIm = handles.currentImage;
    handles.ProjImage = ProjIm(ywin,xwin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(xwin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
    
end
%update handles structure
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


% --- Executes on button press in roiShow.
function roiShow_Callback(hObject, eventdata, handles)
% hObject    handle to roiShow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roiShow
if get(hObject,'Value')
    
    plotROI(handles.Process.ROIpos, handles.StackAxes);
else
    delete(findobj(handles.StackAxes,'LineWidth',1));
end


% --- Executes on selection change in ROIselPop.
function ROIselPop_Callback(hObject, eventdata, handles)
% hObject    handle to ROIselPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ROIselPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ROIselPop
ROI_sel = get(hObject,'Value');
curr_ROI = handles.currentTr(1,5);
if ROI_sel > 1 && curr_ROI ~= ROI_sel - 1
    Trk_tmp = [];
    for i = 1:handles.nTracks
        newTr = handles.outTracks{i};
        Trk_tmp = [Trk_tmp; newTr];
    end
    Trk_tmp = Trk_tmp(Trk_tmp(:,5) == ROI_sel - 1,:);
    if ~isempty(Trk_tmp)
        handles.currentTrIx = Trk_tmp(1,4);
        handles.currentTrIxROI = 1;
        handles.nTracksROI = length(unique(Trk_tmp(:,4)));
    else
        handles.currentTrIxROI = 0;
        handles.nTracksROI = 0;
    end
    
    set(handles.ROITrackCounter,'Visible','on');
    set(handles.ROITrackCounter,'String',['Track ', num2str(handles.currentTrIxROI), '/', num2str(handles.nTracksROI), ' in ROI']);
    
elseif ROI_sel > 1 && curr_ROI == ROI_sel - 1
    Trk_tmp = [];
    for i = 1:handles.nTracks
        newTr = handles.outTracks{i};
        Trk_tmp = [Trk_tmp; newTr];
    end
    Trk_tmp = Trk_tmp(Trk_tmp(:,5) == ROI_sel - 1,:);
    if ~isempty(Trk_tmp)
        
        
        ROITracksNums = unique(Trk_tmp(:,4));
        handles.currentTrIxROI = find(ROITracksNums == handles.currentTrIx);
        handles.nTracksROI = length(unique(Trk_tmp(:,4)));
    else
        handles.currentTrIxROI = 0;
        handles.nTracksROI = 0;
    end
    set(handles.ROITrackCounter,'Visible','on');
    set(handles.ROITrackCounter,'String',['Track ', num2str(handles.currentTrIxROI), '/', num2str(handles.nTracksROI),' in ROI']);
else
    
    
    set(handles.ROITrackCounter,'Visible','off');
end
if handles.currentTrIx > 1
    set(handles.prevTrack,'Enable','on');
end
if handles.currentTrIx < handles.nTracks
    set(handles.nextTrack,'Enable','on');
end
handles.currentTr = handles.outTracks{handles.currentTrIx};

% Set the track counter
set(handles.trackCounter, 'String',...
    ['Track ', num2str(handles.currentTrIx),'/',num2str(handles.nTracks)]);

% Set current Image
handles.currentIx = min(handles.currentTr(:,3));      % current image index

% Determine if we want to display Original or Filtered image
imageType = get(handles.imagePopup,'Value');

if imageType == 1 %Original
    handles.currentImage = ...                            % current image
        handles.Stack(handles.currentIx).data;
elseif imageType == 2 %Filtered
    handles.currentImage = ...
        handles.Process.filterStack(handles.currentIx).data;
end


% Set the image-counter
set(handles.imageCounter, 'String',...
    ['Image ', num2str(handles.currentIx),'/',num2str(handles.nImages)]);

% Set the image-slider
set(handles.ImageSlider,'Value', handles.currentIx);


% Plot first frame of the track in the left plot

axes(handles.StackAxes);
imagesc(handles.currentImage, handles.clims);       % plot image
axis image;

%Plot ROI
if get(handles.roiShow,'Value')
    plotROI(handles.Process.ROIpos,handles.StackAxes);
end

curProj = get(handles.ProjSelPop,'Value');

pos1 = handles.currentTr(1,1:2);
t1 = handles.currentTr(1,3);
tend = handles.currentTr(end,3);
SpExtent = str2double(get(handles.SpaceWindow,'String'));
TiExtent = str2double(get(handles.TimeWindow,'String'));
xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));


% twin = max(1,min(handles.currentIx,t1-TiExtent)):min(handles.nImages,max(handles.currentIx,tend+TiExtent));
% if get(handles.showParticle, 'Value')
%     plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
%     hold on
%     plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
%     hold off
% end
twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    hold on
    plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
    
    hold off
end
if curProj == 1
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,xwin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
    hold on;
    plot([t1,t1],[xwin(1),xwin(end)],'g--');
    plot([tend,tend],[xwin(1),xwin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
    
    
elseif curProj == 2
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
    hold on;
    plot([t1,t1],[ywin(1),ywin(end)],'g--');
    plot([tend,tend],[ywin(1),ywin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
elseif curProj == 3
    ProjIm = handles.currentImage;
    handles.ProjImage = ProjIm(ywin,xwin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(xwin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
    
end


% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ROIselPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROIselPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ProjSelPop.
function ProjSelPop_Callback(hObject, eventdata, handles)
% hObject    handle to ProjSelPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ProjSelPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ProjSelPop

curProj = get(hObject,'Value');

pos1 = handles.currentTr(1,1:2);
t1 = handles.currentTr(1,3);
tend = handles.currentTr(end,3);
SpExtent = str2double(get(handles.SpaceWindow,'String'));
TiExtent = str2double(get(handles.TimeWindow,'String'));
xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));


twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    hold on
    plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
    
    hold off
end
if curProj == 1
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,xwin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
    hold on;
    plot([t1,t1],[xwin(1),xwin(end)],'g--');
    plot([tend,tend],[xwin(1),xwin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
    
    
elseif curProj == 2
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
    hold on;
    plot([t1,t1],[ywin(1),ywin(end)],'g--');
    plot([tend,tend],[ywin(1),ywin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
elseif curProj == 3
    ProjIm = handles.currentImage;
    handles.ProjImage = ProjIm(ywin,xwin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(xwin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
    
end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ProjSelPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ProjSelPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SpaceWindow_Callback(hObject, eventdata, handles)
% hObject    handle to SpaceWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SpaceWindow as text
%        str2double(get(hObject,'String')) returns contents of SpaceWindow as a double
curProj = get(handles.ProjSelPop,'Value');

pos1 = handles.currentTr(1,1:2);
t1 = handles.currentTr(1,3);
tend = handles.currentTr(end,3);
SpExtent = str2double(get(hObject,'String'));
TiExtent = str2double(get(handles.TimeWindow,'String'));
xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));

% twin = max(1,min(handles.currentIx,t1-TiExtent)):min(handles.nImages,max(handles.currentIx,tend+TiExtent));
axes(handles.StackAxes);
% if get(handles.showParticle, 'Value')
%     plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
%     delete(findobj(handles.StackAxes,'Color','w','LineStyle','--'));
%     hold on
%     plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
%     hold off
% end
twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    hold on
    plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
    
    hold off
end
if curProj == 1
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,xwin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
    hold on;
    plot([t1,t1],[xwin(1),xwin(end)],'g--');
    plot([tend,tend],[xwin(1),xwin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
    
    
elseif curProj == 2
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
    hold on;
    plot([t1,t1],[ywin(1),ywin(end)],'g--');
    plot([tend,tend],[ywin(1),ywin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
elseif curProj == 3
    ProjIm = handles.currentImage;
    handles.ProjImage = ProjIm(ywin,xwin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(xwin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
    
end
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function SpaceWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SpaceWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TimeWindow_Callback(hObject, eventdata, handles)
% hObject    handle to TimeWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TimeWindow as text
%        str2double(get(hObject,'String')) returns contents of TimeWindow as a double
curProj = get(handles.ProjSelPop,'Value');

pos1 = handles.currentTr(1,1:2);
t1 = handles.currentTr(1,3);
tend = handles.currentTr(end,3);
SpExtent = str2double(get(handles.SpaceWindow,'String'));
TiExtent = str2double(get(hObject,'String'));
xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));


twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    hold on
    plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
    
    hold off
end
if curProj == 1
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,xwin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
    hold on;
    plot([t1,t1],[xwin(1),xwin(end)],'g--');
    plot([tend,tend],[xwin(1),xwin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
    
    
elseif curProj == 2
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
    hold on;
    plot([t1,t1],[ywin(1),ywin(end)],'g--');
    plot([tend,tend],[ywin(1),ywin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
elseif curProj == 3
    ProjIm = handles.currentImage;
    handles.ProjImage = ProjIm(ywin,xwin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(xwin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
    
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function TimeWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in EndTrackAtT.
function EndTrackAtT_Callback(hObject, eventdata, handles)
% hObject    handle to EndTrackAtT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
posLast = handles.currentTr(end,1:2);
lastT = handles.currentTr(end,3);
curT = handles.currentIx;
refusedMerge = [];
if curT > lastT
    for iFrame = lastT:curT
        FrameInd = find(handles.currentTr(:,3) == iFrame);
        if isempty(FrameInd)
            posLast = handles.currentTr(handles.currentTr(:,3) == iFrame-1,1:2);
            [newPos,intensity,sigma,BG] = singlePeak_FIT(handles.Stack(iFrame).data, posLast);
            
            
            newPos(3) = iFrame;
            newPos(4) = handles.currentTrIx;
            if isfield(handles.Process,'ROIimage')
                newPos = InsideROIcheck2(newPos,handles.Process.ROIimage);
                if isempty(newPos)
                    newPos = handles.currentTr(handles.currentTr(:,3) == iFrame-1,:);
                    newPos(1,3) = iFrame;
                    intensity = newPos(1,3);
                    sigma = 0;
                    BG = 0;
                end
            end
            NewParticle(1) = newPos(1);
            NewParticle(2) = newPos(2);
            
            NewParticle(3) = intensity;
            NewParticle(5) = intensity;
            
            
            NewParticle(6) = newPos(3);
            if size(handles.outParticles,2)> 10;
                
                NewParticle(7) = intensity;
                NewParticle(8) = sigma;
                NewParticle(9) = BG;
                
                NewParticle(10) = newPos(1);
                NewParticle(11) = newPos(2);
                NewParticle(12) = 0;
                NewParticle(13) = newPos(5);
            else
                NewParticle(7) = newPos(5);
            end
            
            
            idxMerge = [];
            
            % check if there's another track with a particle in that
            % position [same frame, within one pixel]
            for iTrack = handles.currentTrIx+1: handles.nTracks
                idxMerge = find(handles.outTracks{1,iTrack}(:,3) == iFrame...
                    & abs(handles.outTracks{1, iTrack}(:,1)-newPos(1))<1 ...
                    & abs(handles.outTracks{1, iTrack}(:,2)-newPos(2))<1);
                
                if ~isempty(idxMerge)
                    break
                end
            end
            
            % if there's already a particle at the selected position,
            %   ask if you want to merge the two tracks
            
            if ~isempty(idxMerge) && isempty(find(refusedMerge == iTrack,1,'first'))
                
                trSegment1 = handles.currentTr;
                trSegment2 = handles.outTracks{1, iTrack}(idxMerge:end,:);
                curProj = get(handles.ProjSelPop,'Value');
                plotTwotrSegments(1,trSegment1, trSegment2, handles.TrackAxes,4-curProj);
                button = questdlg(['Found another track (Track #',...
                    num2str(trSegment2(1,4)),...
                    ') at this position. Do you want to merge the two tracks?'],...
                    'Merge Track','Yes','No','Cancel','Yes');
                plotTwotrSegments(0);
                if strcmp(button,'Cancel')
                    break
                end
                
                if strcmp(button, 'Yes')
                    trSegment2(:,4) = trSegment1(1,4);
                    % Reindex particle identifier
                    for i  = iTrack + 1:handles.nTracks
                        handles.outTracks{i}(:,4) =handles.outTracks{i}(:,4)-1;
                    end
                    % delete track with index iTrack
                    handles.outTracks(iTrack) = [];
                    
                    % update merged track
                    handles.outTracks{handles.currentTrIx} = [trSegment1; trSegment2];
                    handles.currentTr = handles.outTracks{handles.currentTrIx};
                    
                    % update tracks number and track counter
                    
                    handles.nTracks = handles.nTracks-1;
                    set(handles.trackCounter, 'String',...
                        ['Track ', num2str(handles.currentTrIx),'/',num2str(handles.nTracks)]);
                    drawnow
                    if max(handles.currentTr(:,3)) >= curT
                        break
                    end
                else
                    refusedMerge = [refusedMerge; iTrack];
                    PartExist = find(handles.outParticles(:,6) == iFrame...
                        & abs(handles.outParticles(:,1)-newPos(1))<1 ...
                        & abs(handles.outParticles(:,2)-newPos(2))<1);
                    handles.currentTr = [handles.currentTr;newPos];
                    handles.outTracks{handles.currentTrIx} = handles.currentTr;
                    if isempty(PartExist)
                        handles.outParticles = [handles.outParticles; NewParticle];
                    end
                end
                
                %  if there's no other particle at the selected position,
                %  just add the selected position
            else
                PartExist = find(handles.outParticles(:,6) == iFrame...
                    & abs(handles.outParticles(:,1)-newPos(1))<1 ...
                    & abs(handles.outParticles(:,2)-newPos(2))<1);
                
                handles.currentTr = [handles.currentTr;newPos];
                handles.outTracks{handles.currentTrIx} = handles.currentTr;
                if isempty(PartExist)
                    handles.outParticles = [handles.outParticles; NewParticle];
                end
            end
        end
        
    end
elseif curT < lastT
    deleteChoice = questdlg('Tracking data after current point will be deleted. Are you sure you want to continue?','Delete Tracking data','Yes','No','Yes');
    if strcmp(deleteChoice,'Yes')
        idx = find(handles.currentTr(:,3)> curT);
        tindx = handles.currentTr(idx,3);
        % first delete the from Particles
        % update output particles
        if length (handles.outParticles(1,:)) > 10
            PcoordX = 10;
            PcoordY = 11;
        else
            PcoordX = 1;
            PcoordY = 2;
        end
        for iX = 1:length(idx)
            oldPos = handles.currentTr(handles.currentTr(:,3) == tindx(iX),[1 2]);
            pIx = find(handles.outParticles(:,PcoordX) == oldPos(1) & ...
                handles.outParticles(:,PcoordY) == oldPos(2) & ...
                handles.outParticles (:,6) == handles.currentIx);
            handles.outParticles(pIx, :) = [];
            if ~isempty(idx)
                handles.currentTr(handles.currentTr(:,3) == tindx(iX),:) = [];
                handles.outTracks{handles.currentTrIx} = handles.currentTr;
            end
        end
    end
end
set(handles.editTrack, 'Enable', 'on');
curProj = get(handles.ProjSelPop,'Value');

pos1 = handles.currentTr(1,1:2);
t1 = handles.currentTr(1,3);
tend = handles.currentTr(end,3);
SpExtent = str2double(get(handles.SpaceWindow,'String'));
TiExtent = str2double(get(handles.TimeWindow,'String'));
xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));


% twin = max(1,min(handles.currentIx,t1-TiExtent)):min(handles.nImages,max(handles.currentIx,tend+TiExtent));
% if get(handles.showParticle, 'Value')
%     axes(handles.StackAxes);
%     plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
%     hold on
%     plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
%     hold off
% end
twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    hold on
    plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
    
    hold off
end
if curProj == 1
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,xwin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
    hold on;
    plot([t1,t1],[xwin(1),xwin(end)],'g--');
    plot([tend,tend],[xwin(1),xwin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
    
    
elseif curProj == 2
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
    hold on;
    plot([t1,t1],[ywin(1),ywin(end)],'g--');
    plot([tend,tend],[ywin(1),ywin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
elseif curProj == 3
    ProjIm = handles.currentImage;
    handles.ProjImage = ProjIm(ywin,xwin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(xwin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
    
end

guidata(hObject,handles);
outTracks = cell2mat(handles.outTracks');

setappdata(0,'outTracks',outTracks);
setappdata(0,'outParticles',handles.outParticles);


% --- Executes on button press in StartTrackAtT.
function StartTrackAtT_Callback(hObject, eventdata, handles)
% hObject    handle to StartTrackAtT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
posFirst = handles.currentTr(1,1:2);
firstT = handles.currentTr(1,3);
curT = handles.currentIx;
refusedMerge = [];
if curT < firstT
    toAdd = firstT - curT;
    for iFrame = 1:toAdd
        FrameInd = find(handles.currentTr(:,3) == firstT-iFrame);
        if isempty(FrameInd)
            posFirst = handles.currentTr(1,1:2);
            [newPos,intensity,sigma,BG] = singlePeak_FIT(handles.Stack(firstT-iFrame).data, posFirst);
            
            
            
            newPos(3) = firstT-iFrame;
            newPos(4) = handles.currentTrIx;
            if isfield(handles.Process,'ROIimage')
                newPos = InsideROIcheck2(newPos,handles.Process.ROIimage);
            end
            NewParticle(1) = newPos(1);
            NewParticle(2) = newPos(2);
            
            NewParticle(3) = intensity;
            NewParticle(5) = intensity;
            
            
            NewParticle(6) = newPos(3);
            if size(handles.outParticles,2)> 10;
                
                NewParticle(7) = intensity;
                NewParticle(8) = sigma;
                NewParticle(9) = BG;
                
                
                NewParticle(10) = newPos(1);
                NewParticle(11) = newPos(2);
                NewParticle(12) = 0;
                NewParticle(13) = newPos(5);
            else
                NewParticle(7) = newPos(5);
            end
            
            idxMerge = [];
            
            % check if there's another track with a particle in that
            % position [same frame, within one pixel]
            for iTrack = 1: handles.currentTrIx-1
                idxMerge = find(handles.outTracks{1,iTrack}(:,3) == firstT-iFrame...
                    & abs(handles.outTracks{1, iTrack}(:,1)-newPos(1))<1 ...
                    & abs(handles.outTracks{1, iTrack}(:,2)-newPos(2))<1);
                
                if ~isempty(idxMerge)
                    break
                end
            end
            
            % if there's already a particle at the selected position,
            %   ask if you want to merge the two tracks
            
            if ~isempty(idxMerge) && isempty(find(refusedMerge == iTrack,1,'first'))
                
                trSegment1 = handles.outTracks{1, iTrack}(1:idxMerge,:);
                trSegment2 = handles.currentTr;
                trSegment2(:,4) = trSegment1(1,4);
                curProj = get(handles.ProjSelPop,'Value');
                plotTwotrSegments(1,trSegment1, trSegment2, handles.TrackAxes,4-curProj);
                button = questdlg(['Found another track (Track #',...
                    num2str(trSegment1(1,4)),...
                    ') at this position. Do you want to merge the two tracks?'],...
                    'Merge Track','Yes','No','Cancel''Yes');
                plotTwotrSegments(0);
                if strcmp(button,'Cancel')
                    break
                end
                
                if strcmp(button, 'Yes')
                    % Reindex particle identifier
                    for i  = handles.currentTrIx+1:handles.nTracks
                        handles.outTracks{i}(:,4) =handles.outTracks{i}(:,4)-1;
                    end
                    
                    % delete current track
                    handles.outTracks(handles.currentTrIx) =  [];
                    
                    % update merged track
                    handles.outTracks{iTrack} = [trSegment1; trSegment2];
                    handles.currentTr = handles.outTracks{iTrack};
                    handles.currentTrIx = iTrack;
                    
                    % update tracks number and track counter
                    
                    handles.nTracks = handles.nTracks-1;
                    set(handles.trackCounter, 'String',...
                        ['Track ', num2str(handles.currentTrIx),'/',num2str(handles.nTracks)]);
                    drawnow;
                    if min(handles.currentTr(:,3)) <= curT
                        break
                    end
                else
                    PartExist = find(handles.outParticles(:,6) == firstT-iFrame...
                        & abs(handles.outParticles(:,1)-newPos(1))<1 ...
                        & abs(handles.outParticles(:,2)-newPos(2))<1);
                    refusedMerge = [refusedMerge; iTrack];
                    handles.currentTr = [newPos;handles.currentTr];
                    handles.outTracks{handles.currentTrIx} = handles.currentTr;
                    if isempty(PartExist)
                        handles.outParticles = [handles.outParticles; NewParticle];
                    end
                end
                
                %  if there's no other particle at the selected position,
                %  just add the selected position
            else
                PartExist = find(handles.outParticles(:,6) == firstT-iFrame...
                    & abs(handles.outParticles(:,1)-newPos(1))<1 ...
                    & abs(handles.outParticles(:,2)-newPos(2))<1);
                handles.currentTr = [newPos;handles.currentTr];
                handles.outTracks{handles.currentTrIx} = handles.currentTr;
                if isempty(PartExist)
                    handles.outParticles = [NewParticle;handles.outParticles];
                end
            end
        end
        
    end
elseif curT > firstT
    deleteChoice = questdlg('Tracking data after current point will be deleted. Are you sure you want to continue?','Delete Tracking data','Yes','No','Yes');
    if strcmp(deleteChoice,'Yes')
        idx = find(handles.currentTr(:,3)< curT);
        tindx = handles.currentTr(idx,3);
        % first delete the from Particles
        % update output particles
        if length (handles.outParticles(1,:)) > 10
            PcoordX = 10;
            PcoordY = 11;
        else
            PcoordX = 1;
            PcoordY = 2;
        end
        for iX = 1:length(idx)
            oldPos = handles.currentTr(handles.currentTr(:,3) == tindx(iX),[1 2]);
            pIx = find(handles.outParticles(:,PcoordX) == oldPos(1) & ...
                handles.outParticles(:,PcoordY) == oldPos(2) & ...
                handles.outParticles (:,6) == handles.currentIx);
            handles.outParticles(pIx, :) = [];
            if ~isempty(idx)
                handles.currentTr(handles.currentTr(:,3) == tindx(iX),:) = [];
                handles.outTracks{handles.currentTrIx} = handles.currentTr;
            end
        end
%         if ~isempty(idx)
%             handles.currentTr(handles.currentTr(:,3) == tindx(iX),:) = [];
%             handles.outTracks{handles.currentTrIx} = handles.currentTr;
%         end
    end
end
set(handles.editTrack, 'Enable', 'on');
curProj = get(handles.ProjSelPop,'Value');

pos1 = handles.currentTr(1,1:2);
t1 = handles.currentTr(1,3);
tend = handles.currentTr(end,3);
SpExtent = str2double(get(handles.SpaceWindow,'String'));
TiExtent = str2double(get(handles.TimeWindow,'String'));
xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));


% twin = max(1,min(handles.currentIx,t1-TiExtent)):min(handles.nImages,max(handles.currentIx,tend+TiExtent));
% if get(handles.showParticle, 'Value')
%     axes(handles.StackAxes);
%     plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
%     delete(findobj(handles.StackAxes,'Color','w','LineStyle','--'));
%     hold on
%     plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
%     hold off
% end
twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    hold on
    plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');
    
    hold off
end
if curProj == 1
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,xwin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
    hold on;
    plot([t1,t1],[xwin(1),xwin(end)],'g--');
    plot([tend,tend],[xwin(1),xwin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
    
    
elseif curProj == 2
    ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
    handles.ProjImage = ProjIm(:,twin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(twin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
    hold on;
    plot([t1,t1],[ywin(1),ywin(end)],'g--');
    plot([tend,tend],[ywin(1),ywin(end)],'r--');
    plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
    hold off;
    set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
elseif curProj == 3
    ProjIm = handles.currentImage;
    handles.ProjImage = ProjIm(ywin,xwin);
    axes(handles.TrackAxes);
    colormap(gray);                                     % Set colormap
    
    imagesc(xwin,ywin,handles.ProjImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);
    
end

guidata(hObject,handles);
outTracks = cell2mat(handles.outTracks');

setappdata(0,'outTracks',outTracks);
setappdata(0,'outParticles',handles.outParticles);


% --- Executes on button press in KymoSaveButton.
function KymoSaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to KymoSaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'ImSaveDir')
    defaultDir = [handles.ImSaveDir, filesep];
else
    defaultDir = [pwd, filesep];
end

[fname, pname] = uiputfile('*.tif','Save Kymograph',defaultDir);
handles.ImSaveDir = pname;

if fname ~= 0
    imwrite(handles.ProjImage,[pname, filesep, fname]);
end
guidata(hObject,handles);


% --- Executes on button press in FirstFrameSaveButton.
function FirstFrameSaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to FirstFrameSaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'ImSaveDir')
    defaultDir = [handles.ImSaveDir, filesep];
else
    defaultDir = [pwd, filesep];
end

[fname, pname] = uiputfile('*.tif','Save First Frame',defaultDir);
handles.ImSaveDir = pname;
imageType = get(handles.imagePopup,'Value');

if imageType == 1 %Original
    FirstImage = ...                            % current image
        handles.Stack(1).data;
elseif imageType == 2 %Filtered
    FirstImage = ...
        uint16(handles.Process.filterStack(1).data);
end

if fname ~= 0
    imwrite(FirstImage,[pname, filesep, fname]);
end
guidata(hObject,handles);


% --- Executes on button press in MoveToLastPos.
function MoveToLastPos_Callback(hObject, eventdata, handles)
% hObject    handle to MoveToLastPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idx = find(handles.currentTr(:,3) == handles.currentIx);

idx_Prev = find(handles.currentTr(:,3) == handles.currentIx - 1);
idx_Next = find(handles.currentTr(:,3) == handles.currentIx + 1);

if ~isempty(idx_Prev) && ~isempty(idx_Next)
    DirectChoice = questdlg('Do you want to use the previous or the next point from the track?','Choose Direction','Previous','Next','Previous');
elseif ~isempty(idx_Prev) && isempty(idx_Next)
    DirectChoice = 'Previous';
elseif isempty(idx_Prev) && ~isempty(idx_Next)
    DirectChoice = 'Next';
else
    DirectChoice = '';
end

if strcmp(DirectChoice,'Next')
    pos = handles.currentTr(idx_Next,1:2);
else
    pos = handles.currentTr(idx_Prev,1:2);
end
if ~strcmp(DirectChoice,'');
    FitGauss = get(handles.gaussianFit, 'Value');
    if isempty(idx)
        newPos(1:2) = pos;
        
        
        
        if FitGauss % refine position with Gaussian fitting
            
            [newPos,intensity,sigma,BG] = singlePeak_FIT(handles.Stack(handles.currentIx).data, newPos);
            
        
            NewParticle(3) = intensity;
            NewParticle(5) = intensity;

        else
            NewParticle(3) = handles.currentImage(round(newPos(2)),round(newPos(1)));
            NewParticle(5) = handles.currentImage(round(newPos(2)),round(newPos(1)));
        end
        newPos(3) = handles.currentIx;
        newPos(4) = handles.currentTrIx;
        if isfield(handles.Process,'ROIimage')
            newPos = InsideROIcheck2(newPos,handles.Process.ROIimage);
        end
        NewParticle(1) = newPos(1);
        NewParticle(2) = newPos(2);
        NewParticle(6) = newPos(3);
        if size(handles.outParticles,2)> 10;
            if FitGauss
                NewParticle(7) = intensity;
                NewParticle(8) = sigma;
                NewParticle(9) = BG;
            end
            NewParticle(10) = newPos(1);
            NewParticle(11) = newPos(2);
            NewParticle(12) = 0;
            NewParticle(13) = newPos(5);
        else
            NewParticle(7) = newPos(5);
        end



        % check if the particle is added at the beginning of the track
        if handles.currentIx == min(handles.currentTr(:,3)) - 1;
            idxMerge = [];

            % check if there's another track with a particle in that
            % position [same frame, within one pixel]
            for iTrack = 1:handles.currentTrIx-1
                idxMerge = find(handles.outTracks{1,iTrack}(:,3) == handles.currentIx...
                    & abs(handles.outTracks{1, iTrack}(:,1)-newPos(1))<1 ...
                    & abs(handles.outTracks{1, iTrack}(:,2)-newPos(2))<1);

                if ~isempty(idxMerge)
                    break
                end
            end

            % if there's already a particle at the selected position,
            %   ask if you want to merge the two tracks

            if ~isempty(idxMerge)

                trSegment1 = handles.outTracks{1, iTrack}(1:idxMerge,:);
                trSegment2 = handles.currentTr;
                trSegment2(:,4) = trSegment1(1,4);

                plotTwotrSegments(1,trSegment1, trSegment2, handles.TrackAxes,(4-get(handles.ProjSelPop,'Value')));
                button = questdlg(['Found another track (Track #',...
                    num2str(trSegment1(1,4)),...
                    ') at this position. Do you want to merge the two tracks?'],...
                    'Merge Track','Yes','No','Yes');
                plotTwotrSegments(0);

                if strcmp(button, 'Yes')
                    % Reindex particle identifier
                    for i  = handles.currentTrIx+1:handles.nTracks
                        handles.outTracks{i}(:,4) =handles.outTracks{i}(:,4)-1;
                    end

                    % delete current track
                    handles.outTracks(handles.currentTrIx) =  [];

                    % update merged track
                    handles.outTracks{iTrack} = [trSegment1; trSegment2];
                    handles.currentTr = handles.outTracks{handles.currentTrIx};

                    % update tracks number and track counter

                    handles.nTracks = handles.nTracks-1;
                    set(handles.trackCounter, 'String',...
                        ['Track ', num2str(handles.currentTrIx),'/',num2str(handles.nTracks)]);
                else
                    PartExist = find(handles.outParticles(:,6) == handles.currentIx...
                        & abs(handles.outParticles(:,1)-newPos(1))<1 ...
                        & abs(handles.outParticles(:,2)-newPos(2))<1);
                    handles.currentTr = [newPos;handles.currentTr];
                    handles.outTracks{handles.currentTrIx} = handles.currentTr;
                    if isempty(PartExist)
                        handles.outParticles = [handles.outParticles; NewParticle];
                    end

                end

            else % if no merging is possible just add the particle
                PartExist = find(handles.outParticles(:,6) == handles.currentIx...
                    & abs(handles.outParticles(:,1)-newPos(1))<1 ...
                    & abs(handles.outParticles(:,2)-newPos(2))<1);
                handles.currentTr = [newPos;handles.currentTr];
                handles.outTracks{handles.currentTrIx} = handles.currentTr;
                if isempty(PartExist)
                    handles.outParticles = [handles.outParticles; NewParticle];
                end
            end

            % check if the particle is added at the beginning of the track

        elseif handles.currentIx == max(handles.currentTr(:,3)) + 1;

            idxMerge = [];

            % check if there's another track with a particle in that
            % position [same frame, within one pixel]
            for iTrack = handles.currentTrIx+1: handles.nTracks
                idxMerge = find(handles.outTracks{1,iTrack}(:,3) == handles.currentIx...
                    & abs(handles.outTracks{1, iTrack}(:,1)-newPos(1))<1 ...
                    & abs(handles.outTracks{1, iTrack}(:,2)-newPos(2))<1);

                if ~isempty(idxMerge)
                    break
                end
            end

            % if there's already a particle at the selected position,
            %   ask if you want to merge the two tracks

            if ~isempty(idxMerge)

                trSegment1 = handles.currentTr;
                trSegment2 = handles.outTracks{1, iTrack}(idxMerge:end,:);

                plotTwotrSegments(1,trSegment1, trSegment2, handles.TrackAxes,(4-get(handles.ProjSelPop,'Value')));
                button = questdlg(['Found another track (Track #',...
                    num2str(trSegment2(1,4)),...
                    ') at this position. Do you want to merge the two tracks?'],...
                    'Merge Track','Yes','No','Yes');
                plotTwotrSegments(0);

                if strcmp(button, 'Yes')
                    trSegment2(:,4) = trSegment1(1,4);
                    % Reindex particle identifier
                    for i  = iTrack + 1:handles.nTracks
                        handles.outTracks{i}(:,4) =handles.outTracks{i}(:,4)-1;
                    end
                    % delete track with index iTrack
                    handles.outTracks(iTrack) = [];

                    % update merged track
                    handles.outTracks{handles.currentTrIx} = [trSegment1; trSegment2];
                    handles.currentTr = handles.outTracks{handles.currentTrIx};

                    % update tracks number and track counter

                    handles.nTracks = handles.nTracks-1;
                    set(handles.trackCounter, 'String',...
                        ['Track ', num2str(handles.currentTrIx),'/',num2str(handles.nTracks)]);

                end

                %  if there's no other particle at the selected position,
                %  just add the selected position
            else
                handles.currentTr = [handles.currentTr;newPos];
                handles.outTracks{handles.currentTrIx} = handles.currentTr;
                handles.outParticles = [handles.outParticles; NewParticle];
            end

        else
            msgbox('You can only add particles to frames close to the current track!');
        end

        %Plot ROI
        if get(handles.roiShow,'Value')
            plotROI(handles.Process.ROIpos,handles.StackAxes);
        end
        % Plot the current Track
        plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
        curProj = get(handles.ProjSelPop,'Value');

        pos1 = handles.currentTr(1,1:2);
        t1 = handles.currentTr(1,3);
        tend = handles.currentTr(end,3);
        SpExtent = str2double(get(handles.SpaceWindow,'String'));
        TiExtent = str2double(get(handles.TimeWindow,'String'));
        xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
        ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));


        twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
        if get(handles.showParticle, 'Value')
            plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
            hold on
            plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');

            hold off
        end
        if curProj == 1
            ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
            handles.ProjImage = ProjIm(:,twin);

            axes(handles.TrackAxes);
            colormap(gray);                                     % Set colormap

            imagesc(twin,xwin,handles.ProjImage, handles.clims);
            axis image;

            plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
            hold on;
            plot([t1,t1],[xwin(1),xwin(end)],'g--');
            plot([tend,tend],[xwin(1),xwin(end)],'r--');
            plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
            hold off;
            set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);


        elseif curProj == 2
            ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
            handles.ProjImage = ProjIm(:,twin);
            axes(handles.TrackAxes);
            colormap(gray);                                     % Set colormap

            imagesc(twin,ywin,handles.ProjImage, handles.clims);
            axis image;

            plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
            hold on;
            plot([t1,t1],[ywin(1),ywin(end)],'g--');
            plot([tend,tend],[ywin(1),ywin(end)],'r--');
            plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
            hold off;
            set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
        elseif curProj == 3
            ProjIm = handles.currentImage;
            handles.ProjImage = ProjIm(ywin,xwin);
            axes(handles.TrackAxes);
            colormap(gray);                                     % Set colormap

            imagesc(xwin,ywin,handles.ProjImage, handles.clims);
            axis image;

            plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);

        end

        set(handles.editTrack,'Enable','on');

    else %move particle
        newPos = pos;
        if FitGauss % refine position with Gaussian fitting

            newPos = singlePeak_FIT(handles.Stack(handles.currentIx).data, newPos);

        end

        % update output track
        idx = find(handles.currentTr(:,3) == handles.currentIx);
        oldPos = handles.currentTr(idx, [1 2]);
        handles.currentTr(idx,[1 2]) = newPos;
        handles.outTracks{handles.currentTrIx} = handles.currentTr;

        % update output particles
        if length (handles.outParticles(1,:)) > 10
            PcoordX = 10;
            PcoordY = 11;
            ROIindx = 13;
        else
            PcoordX = 1;
            PcoordY = 2;
            ROIindx = 7;
        end
        pIx = find(handles.outParticles(:,PcoordX) == oldPos(1) & ...
            handles.outParticles(:,PcoordY) == oldPos(2) & ...
            handles.outParticles (:,6) == handles.currentIx);

        handles.outParticles(pIx,PcoordX) = newPos(1);
        handles.outParticles(pIx,PcoordY) = newPos(2);
        if isfield(handles.Process,'ROIimage')
            newPos = InsideROIcheck2([newPos handles.currentIx handles.currentTrIx],handles.Process.ROIimage);
            handles.outParticles(pIx,ROIindx) = newPos(5);
        end

        %Plot ROI
        if get(handles.roiShow,'Value')
            plotROI(handles.Process.ROIpos,handles.StackAxes);
        end

        % plot new track
        plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
        curProj = get(handles.ProjSelPop,'Value');

        pos1 = handles.currentTr(1,1:2);
        t1 = handles.currentTr(1,3);
        tend = handles.currentTr(end,3);
        SpExtent = str2double(get(handles.SpaceWindow,'String'));
        TiExtent = str2double(get(handles.TimeWindow,'String'));
        xwin = max(1,round(pos1(1)-(SpExtent/2))):min(size(handles.curImTimesWin,2),round(pos1(1)+(SpExtent/2)));
        ywin = max(1,round(pos1(2)-(SpExtent/2))):min(size(handles.curImTimesWin,1),round(pos1(2)+(SpExtent/2)));


        twin = max(1,min(handles.currentIx-round(TiExtent/2),handles.nImages-TiExtent)):min(handles.nImages,max(handles.currentIx+round(TiExtent/2),1+TiExtent));
        if get(handles.showParticle, 'Value')
            plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
            hold on
            plot([xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],[ywin(1) ywin(end) ywin(end) ywin(1) ywin(1)],'w--');

            hold off
        end
        if curProj == 1
            ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[2,3,1]),[],3);
            handles.ProjImage = ProjIm(:,twin);

            axes(handles.TrackAxes);
            colormap(gray);                                     % Set colormap

            imagesc(twin,xwin,handles.ProjImage, handles.clims);
            axis image;

            plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,2);
            hold on;
            plot([t1,t1],[xwin(1),xwin(end)],'g--');
            plot([tend,tend],[xwin(1),xwin(end)],'r--');
            plot([handles.currentIx, handles.currentIx],[xwin(1),xwin(end)],'w');
            hold off;
            set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);


        elseif curProj == 2
            ProjIm = max(permute(handles.curImTimesWin(ywin,xwin,:),[1,3,2]),[],3);
            handles.ProjImage = ProjIm(:,twin);
            axes(handles.TrackAxes);
            colormap(gray);                                     % Set colormap

            imagesc(twin,ywin,handles.ProjImage, handles.clims);
            axis image;

            plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,3);
            hold on;
            plot([t1,t1],[ywin(1),ywin(end)],'g--');
            plot([tend,tend],[ywin(1),ywin(end)],'r--');
            plot([handles.currentIx, handles.currentIx],[ywin(1),ywin(end)],'w');
            hold off;
            set(handles.TrackAxes,'Xlim',[twin(1),twin(end)]);
        elseif curProj == 3
            ProjIm = handles.currentImage;
            handles.ProjImage = ProjIm(ywin,xwin);
            axes(handles.TrackAxes);
            colormap(gray);                                     % Set colormap

            imagesc(xwin,ywin,handles.ProjImage, handles.clims);
            axis image;

            plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes,1);

        end

        % Update handles structure


    end
    guidata(hObject, handles);
    outTracks = cell2mat(handles.outTracks');

    setappdata(0,'outTracks',outTracks);
    setappdata(0,'outParticles',handles.outParticles);
end
