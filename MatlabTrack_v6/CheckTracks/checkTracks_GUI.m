function varargout = checkTracks_GUI(varargin)
% CHECKTRACKS_GUI M-file for checkTracks_GUI.fig
%      CHECKTRACKS_GUI, by itself, creates a new CHECKTRACKS_GUI or raises the existing
%      singleton*.
%
%      H = CHECKTRACKS_GUI returns the handle to a new CHECKTRACKS_GUI or the handle to
%      the existing singleton*.
%
%      CHECKTRACKS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHECKTRACKS_GUI.M with the given input arguments.
%
%      CHECKTRACKS_GUI('Property','Value',...) creates a new CHECKTRACKS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before checkTracks_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to checkTracks_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help checkTracks_GUI

% Last Modified by GUIDE v2.5 24-Jul-2015 10:35:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @checkTracks_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @checkTracks_GUI_OutputFcn, ...
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




% --- Executes just before checkTracks_GUI is made visible.
% ------------------------------------------------------------------------
function checkTracks_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to checkTracks_GUI (see VARARGIN)

% Choose default command line output for checkTracks_GUI


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
MaxBrightnessStack = 2*max(MaxBrightness); %give a 10% buffer above the maximum
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
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
end


% Plot First frame + track overlay in the right plot

axes(handles.TrackAxes);
colormap(gray);                                     % Set colormap
handles.clims = [min(min(handles.currentImage)) ... % Set color scale
max(max(handles.currentImage))];

%handles.clims = [2000, 7000];
imagesc(handles.currentImage, handles.clims);
axis image;

plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes checkTracks_GUI wait for user response (see UIRESUME)
 uiwait(handles.figure1);
% ------------------------------------------------------------------------



% --- Outputs from this function are returned to the command line.
function varargout = checkTracks_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get output from handles structure
varargout{1} = handles.out2{1};
varargout{2} = handles.out2{2};

% The figure can be deleted now
delete(handles.figure1);





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
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
end


axes(handles.TrackAxes);
imagesc(handles.currentImage, handles.clims);
axis image;

plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes);


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

if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
end


% Plot First frame + track overlay in the right plot

axes(handles.TrackAxes);
imagesc(handles.currentImage, handles.clims);
axis image;

plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes);


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

if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
end


% Plot First frame + track overlay in the right plot

axes(handles.TrackAxes);
imagesc(handles.currentImage, handles.clims);
axis image;

plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes);


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
    [newPos(1) newPos(2)] = ginput(1);
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

   newPos = singlePeak_FIT(handles.currentImage, newPos);
    
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
plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes);

% Update handles structure
guidata(hObject, handles);

 




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
    
    [newPos(1) newPos(2)] = ginput(1); 
    
    
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
        
        newPos = singlePeak_FIT(handles.currentImage, newPos);
    
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
       
          plotTwotrSegments(1,trSegment1, trSegment2, handles.TrackAxes);
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
               
            end
            
        else % if no merging is possible just add the particle
            handles.currentTr = [newPos;handles.currentTr];
            handles.outTracks{handles.currentTrIx} = handles.currentTr;
            handles.outParticles = [handles.outParticles; NewParticle];
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
    
          plotTwotrSegments(1,trSegment1, trSegment2, handles.TrackAxes);
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
     plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes);

     
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
        plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes);
    
    else  % if the segment to cut is in the middle of the track Split the particle in two.
       
       trSegment1 = handles.currentTr(1:idx-1,:);
       trSegment2 = handles.currentTr(idx+1:end,:);
       trSegment2(:,4) = trSegment2(:,4) + 1;
       
       plotTwotrSegments(1,trSegment1, trSegment2, handles.TrackAxes);
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
                plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes);

               
           
      
            
               
           

    end 
    % Update handles structure
    guidata(hObject, handles);
end





    


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
                pIx(indx) = find(handles.outParticles(:,PcoordX) == handles.currentTr(i,1) & ...
                    handles.outParticles(:,PcoordY) == handles.currentTr(i,2) & ...
                    handles.outParticles (:,6) == handles.currentTr(i,3));
                indx = indx+1;
            end
                
        end
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
                 if get(handles.showParticle, 'Value')
                     plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
                 end
                 
                 plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes);

                % Update handles structure
                guidata(hObject, handles);
                            
 end


% -------------------------------------------------------------------------




% --- Executes on button press in done.
% -------------------------------------------------------------------------
function done_Callback(hObject, eventdata, handles)

% merge the tracks together
outTracks = cell2mat(handles.outTracks');

% outputtracks
handles.out2{1} = outTracks;
handles.out2{2} = handles.outParticles;

%update handles
guidata(hObject, handles);

% Resume UI
uiresume(handles.figure1);






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
else
    delete(findobj(handles.StackAxes,'Color','r'));
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
    if get(handles.showParticle, 'Value')
        plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
    end
    
    
    % Plot First frame + track overlay in the right plot
    
    axes(handles.TrackAxes);
    imagesc(handles.currentImage, handles.clims);
    axis image;
    
    plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes);
    
    
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

if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
end


% Plot First frame + track overlay in the right plot

axes(handles.TrackAxes);
colormap(gray);                                     % Set colormap

imagesc(handles.currentImage, handles.clims);
axis image;

plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes);


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

% merge the tracks together
outTracks = cell2mat(handles.outTracks');

% outputtracks
handles.out2{1} = outTracks;
handles.out2{2} = handles.outParticles;


%update handles
guidata(hObject, handles);

% Resume UI
uiresume(handles.figure1);


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
          
           
elseif imageType == 2 %Filtered
    handles.currentImage = ...
        handles.Process.filterStack(handles.currentIx).data;
    for i = 1:handles.nImages
        MaxBrightness(i) = max(handles.Process.filterStack(i).data(:));
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
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
end


axes(handles.TrackAxes);
imagesc(handles.currentImage, handles.clims);
axis image;

plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes);

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
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
end


axes(handles.TrackAxes);
imagesc(handles.currentImage, handles.clims);
axis image;

plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes);

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
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
end


axes(handles.TrackAxes);
imagesc(handles.currentImage, handles.clims);
axis image;

plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes);

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
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
end


axes(handles.TrackAxes);
imagesc(handles.currentImage, handles.clims);
axis image;

plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes);

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
BlackLevel = str2double(get(handles.WhiteValEdit,'String'));

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
if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
end


axes(handles.TrackAxes);
imagesc(handles.currentImage, handles.clims);
axis image;

plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes);

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
    handles.currentTrIx = Trk_tmp(1,4);
    
        
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

if get(handles.showParticle, 'Value')
    plotSingleParticle(handles.currentTr,handles.currentIx,handles.StackAxes);
end


% Plot First frame + track overlay in the right plot

axes(handles.TrackAxes);
imagesc(handles.currentImage, handles.clims);
axis image;

plotSingleTrack (handles.currentTr, handles.currentIx, handles.TrackAxes);


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
