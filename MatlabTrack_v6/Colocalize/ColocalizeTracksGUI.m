function varargout = ColocalizeTracksGUI(varargin)
% COLOCALIZETRACKSGUI MATLAB code for ColocalizeTracksGUI.fig
%      COLOCALIZETRACKSGUI, by itself, creates a new COLOCALIZETRACKSGUI or raises the existing
%      singleton*.
%
%      H = COLOCALIZETRACKSGUI returns the handle to a new COLOCALIZETRACKSGUI or the handle to
%      the existing singleton*.
%
%      COLOCALIZETRACKSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COLOCALIZETRACKSGUI.M with the given input arguments.
%
%      COLOCALIZETRACKSGUI('Property','Value',...) creates a new COLOCALIZETRACKSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ColocalizeTracksGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ColocalizeTracksGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ColocalizeTracksGUI

% Last Modified by GUIDE v2.5 29-Jan-2012 16:17:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ColocalizeTracksGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ColocalizeTracksGUI_OutputFcn, ...
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


% --- Executes just before ColocalizeTracksGUI is made visible.
function ColocalizeTracksGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ColocalizeTracksGUI (see VARARGIN)

% Choose default command line output for ColocalizeTracksGUI
handles.output = hObject;


% Get Input

handles.Tracks = varargin{1};
handles.pixelSize = varargin{2};
handles.frameTime = varargin{3};
handles.pathName = varargin{4};

% Load reference image/stack and calculate the average
[handles.fileName,handles.pathName]  = ...
    uigetfile('*.tif','Select reference stack for read out', handles.pathName);

[imageStack, nImages] = TIFread([handles.pathName,handles.fileName]);

handles.SumImage(1).data = double(imageStack(1).data);


for imageIx =  2:nImages;
    handles.SumImage(1).data = handles.SumImage(1).data + ...
        double(imageStack(imageIx).data);
end;
 handles.SumImage(1).data =  handles.SumImage(1).data / nImages;


handles.SumImage(2).data = 0;
handles.projlims = [min(min(handles.SumImage(1).data)) ...
    max(max(handles.SumImage(1).data))];


% Display Images and tracks
axes(handles.PolII_image);
imagesc(handles.SumImage(1).data, handles.projlims);

axis image;
colormap(gray);

% plot tracks;
hold on;

delete(findobj(gca,'Color','y')); 
for i = 1:1:max(handles.Tracks(:,4))
    ix = find(handles.Tracks(:,4)==i);
    plot(handles.Tracks(ix,1),handles.Tracks(ix,2),'g', 'LineWidth',1);
 
end;
hold off;

% Enable the plot tracks checkbox
set(handles.showTracks, 'Value', 1);
set(handles.showTracks, 'Enable', 'on');

% Disable the show particles checkbox
set(handles.showParticles, 'Value', 0);
set(handles.showParticles, 'Enable', 'off');

% Disable the residence time button
set(handles.Res_vs_d,'Enable','off');

% Disable the plot histogram button
set(handles.plotHistogram,'Enable','off');

% Disable the find particles button
set(handles.FindParticles,'Enable','off');



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ColocalizeTracksGUI wait for user response (see UIRESUME)



% --- Outputs from this function are returned to the command line.
function varargout = ColocalizeTracksGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in showTracks.
function showTracks_Callback(hObject, eventdata, handles)
% hObject    handle to showTracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showTracks


isOn = get(hObject, 'Value');
if isOn;
    
   hold on;
    for i = 1:1:max(handles.Tracks(:,4))
        ix = find(handles.Tracks(:,4)==i);
        plot(handles.Tracks(ix,1),handles.Tracks(ix,2),'g', 'LineWidth',1);
 
    end;
    
    hold off;
    
    
else
   delete(findobj(gca,'Color','g')); 
   delete(findobj(gca,'Color','y')); 
end;



% --- Executes on button press in showParticles.
function showParticles_Callback(hObject, eventdata, handles)
% hObject    handle to showParticles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showParticles

isOn = get(hObject, 'Value');
if isOn;

   plotParticle(handles.Centroids, 1, 0);
 
else
    delete(findobj(gca,'Color','r'));
end;


% --- Executes on button press in Res_vs_d.
function Res_vs_d_Callback(hObject, eventdata, handles)
% hObject    handle to Res_vs_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tracks_small = TrStartEnd(handles.Tracks);
distance = distanceFromPol(Tracks_small,handles.Centroids);
Tracks_small(:,7) = distance;

Tracks_small(:,5) = Tracks_small(:,5)*handles.frameTime;
Tracks_small(:,7) = Tracks_small(:,7)*handles.pixelSize;

figure;
plot(Tracks_small(:,7),Tracks_small(:,5),'ok');
xlabel('Distance from Pol [\mum]');
ylabel('Residence Time [s]');

num2clip(Tracks_small);

handles.Tracks_small = Tracks_small;
set(handles.plotHistogram,'Enable','on');

% Update handles structure
guidata(hObject, handles);



function MaxDistance_Callback(hObject, eventdata, handles)
% hObject    handle to MaxDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxDistance as text
%        str2double(get(hObject,'String')) returns contents of MaxDistance as a double


% --- Executes during object creation, after setting all properties.
function MaxDistance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotHistogram.
function plotHistogram_Callback(hObject, eventdata, handles)
% hObject    handle to plotHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(findobj(gca,'Color','y')); 
maxDistance = str2double(get(handles.MaxDistance,'String'));
handles.outPar(5) = maxDistance;
HistResT(:,1) = (1:1:400)*handles.frameTime;
HistResT(:,2) = hist(handles.Tracks_small(:,5),HistResT(:,1));
Test = handles.Tracks_small(:,7) < maxDistance;

hold on;
    for i = 1:1:max(handles.Tracks(:,4))
        if Test(i)
        ix = find(handles.Tracks(:,4)==i);
        plot(handles.Tracks(ix,1),handles.Tracks(ix,2),'y', 'LineWidth',1);
        end
    end;
hold off;

HistResT(:,3) = hist(handles.Tracks_small(Test, 5),HistResT(:,1));

figure;
hold on;
plot(HistResT(:,1),HistResT(:,2),'ok');
plot(HistResT(:,1), HistResT(:,3),'or');
hold off;

xlabel('Residence Time [s]');
ylabel('Counts');

num2clip(HistResT);




function Low_BP_Callback(hObject, eventdata, handles)
% hObject    handle to Low_BP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Low_BP as text
%        str2double(get(hObject,'String')) returns contents of Low_BP as a double


% --- Executes during object creation, after setting all properties.
function Low_BP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Low_BP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function High_BP_Callback(hObject, eventdata, handles)
% hObject    handle to High_BP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of High_BP as text
%        str2double(get(hObject,'String')) returns contents of High_BP as a double


% --- Executes during object creation, after setting all properties.
function High_BP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to High_BP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FilterStack.
function FilterStack_Callback(hObject, eventdata, handles)
% hObject    handle to FilterStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

loBP = str2double(get(handles.Low_BP,'String'));
hiBP = str2double(get(handles.High_BP,'String'));

handles.filterImage(1).data = ...
    bpass(handles.SumImage(1).data, loBP, hiBP);
handles.filterImage(2).data = 0;
handles.outPar(1) = loBP;
handles.outPar(2) = hiBP;

% Enable the find particles button
set(handles.FindParticles,'Enable','on');

set(handles.Res_vs_d,'Enable','off');

set(handles.plotHistogram,'Enable','off');


% Update handles structure
guidata(hObject, handles);



function Threshold_Callback(hObject, eventdata, handles)
% hObject    handle to Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Threshold as text
%        str2double(get(hObject,'String')) returns contents of Threshold as a double


% --- Executes during object creation, after setting all properties.
function Threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ws_Callback(hObject, eventdata, handles)
% hObject    handle to ws (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ws as text
%        str2double(get(hObject,'String')) returns contents of ws as a double


% --- Executes during object creation, after setting all properties.
function ws_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ws (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FindParticles.
function FindParticles_Callback(hObject, eventdata, handles)
% hObject    handle to FindParticles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Threshold = str2double(get(handles.Threshold,'String'));
hiBP = str2double(get(handles.High_BP,'String'));
windowSz = str2double(get(handles.ws,'String'));

handles.outPar(3) = Threshold;
handles.outPar(4) = windowSz;

% Calculate Particles centroids
Centroid = findParticles(handles.filterImage, Threshold, hiBP,windowSz);
Centroid(:,7) = ones(size(Centroid,1),1);
handles.Centroids = Centroid;

plotParticle(Centroid,1,0);
set(handles.showParticles, 'Value', 1);
set(handles.showParticles, 'Enable', 'on');


% Enable the residence time button
set(handles.Res_vs_d,'Enable','on');

set(handles.plotHistogram,'Enable','off');

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in ExportPar.
function ExportPar_Callback(hObject, eventdata, handles)
% hObject    handle to ExportPar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

num2clip(handles.outPar);

% --- Executes on button press in SaveImage.
function SaveImage_Callback(hObject, eventdata, handles)
% hObject    handle to SaveImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

oldax = gca;
fig2 = figure;
newax = copyobj(oldax(1),fig2); 
colormap(gray);
axis off;
[FileName,PathName] = uiputfile('.tif') ;
saveas(fig2,[PathName,FileName]);
close(fig2);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

delete(hObject);
