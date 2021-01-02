function varargout = ChannelAlign_GUI(varargin)
% CHANNELALIGN_GUI MATLAB code for ChannelAlign_GUI.fig
%      CHANNELALIGN_GUI, by itself, creates a new CHANNELALIGN_GUI or raises the existing
%      singleton*.
%
%      H = CHANNELALIGN_GUI returns the handle to a new CHANNELALIGN_GUI or the handle to
%      the existing singleton*.
%
%      CHANNELALIGN_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHANNELALIGN_GUI.M with the given input arguments.
%
%      CHANNELALIGN_GUI('Property','Value',...) creates a new CHANNELALIGN_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ChannelAlign_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ChannelAlign_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ChannelAlign_GUI

% Last Modified by GUIDE v2.5 11-May-2017 12:54:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ChannelAlign_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ChannelAlign_GUI_OutputFcn, ...
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


% --- Executes just before ChannelAlign_GUI is made visible.
function ChannelAlign_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ChannelAlign_GUI (see VARARGIN)

% Choose default command line output for ChannelAlign_GUI
handles.output = hObject;

handles.Images = [];
handles.Particles = [];
handles.Tracks = [];
handles.ImagesAlign.Image1 = [];
handles.ImagesAlign.Image2 = [];
handles.ParticlesAlign.Particles1 = [];
handles.ParticlesAlign.Particles2 = [];
handles.TracksAlign.Tracks1 = [];
handles.TracksAlign.Tracks2 = [];

%get info to set up the default calibration file name
MatPath = path;

MatTrackInPath = strfind(MatPath,'MatlabTrack_');
semicol_char = pathsep;
semicol_loc = strfind(MatPath,semicol_char);

if semicol_loc(1) > MatTrackInPath(1)
    MatTrackPath = MatPath(1:semicol_loc(1)-1);
else
    startInd = find(semicol_loc < MatTrackInPath(1),1,'last');

    endInd = find(semicol_loc > MatTrackInPath(1),1,'first');
    MatTrackPath = MatPath(semicol_loc(startInd)+1:semicol_loc(endInd)-1);
end

CalibPath = [MatTrackPath,filesep,'TwoColorSMcoloc', filesep, 'Calibration Files'];
curFold = pwd;
cd(CalibPath);
CalibFiles = dir('*.mat');

CalibFileNames = cell(size(CalibFiles));



for i = 1:size(CalibFiles)
    CalibFileNames{i,:} = CalibFiles(i).name;
    dates(i) = datenum(CalibFiles(i).date);
end

%Detemine most recent file
mostRecent = find(dates == max(dates));

set(handles.RegPop,'String',CalibFileNames);
set(handles.RegPop,'Value',mostRecent);

%Load in this file
IN = load(CalibFiles(mostRecent).name);

cd(curFold);

handles.CalibData.RegMat = IN.RegMat;
handles.CalibData.wavelengths = IN.wavelengths;
handles.CalibData.CalibPath = CalibPath;
handles.CalibData.Files = CalibFileNames;

%If the GUI was initiated from the 2-color SMC GUI, then there will be
%images, particles and tracks to load in
if ~isempty(varargin)
    handles.Images = varargin{1};
    handles.Particles = varargin{2};
    handles.Tracks = varargin{3};
   
    handles.CurrentInd = 1;
    %display the merged image
    clims1 = [min(handles.Images.Image1(1).data(:)),max(handles.Images.Image1(1).data(:))];
    clims2 =[min(handles.Images.Image2(1).data(:)),max(handles.Images.Image2(1).data(:))];
    
    handles.clims1 = clims1;
    handles.clims2 = clims2;
    im = merge2colorData(handles.Images.Image1(1).data,handles.Images.Image2(1).data,[1 2],clims1,clims2);
    set(handles.AlignButton,'Enable','on');
    
    maxIm = min(size(handles.Images.Image1,2),size(handles.Images.Image2,2));
    set(handles.ImageSlider,'Enable','on');
    set(handles.ImageSlider,'Min',1,'Max',maxIm,'SliderStep',[1/(maxIm - 1), 1/(maxIm - 1)]);
    set(handles.ImageSlider,'Value',1);
    axes(handles.OrigAxes);
    image(im);
    
end

% Update handles structure
guidata(hObject, handles);

setappdata(0,'imageStack1',handles.Images.Image1);
setappdata(0,'imageStack2',handles.Images.Image2);
setappdata(0,'Particles1',handles.Particles.Particles1);
setappdata(0,'Particles2',handles.Particles.Particles2);
setappdata(0,'Tracks1',handles.Tracks.Tracks1);
setappdata(0,'Tracks2',handles.Tracks.Tracks1);


% --- Outputs from this function are returned to the command line.
function varargout = ChannelAlign_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
varargout{1} = handles.output;

% Get default command line output from handles structure



% --- Executes on selection change in RegPop.
function RegPop_Callback(hObject, eventdata, handles)
% hObject    handle to RegPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns RegPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RegPop

%Get the string selected

AllCalibFiles = get(hObject,'String');
SelCalib = get(hObject,'Value');

IN = load([handles.CalibData.CalibPath,filesep,AllCalibFiles{SelCalib,:}]);

handles.CalibData.RegMat = IN.RegMat;
handles.CalibData.wavelengths = IN.wavelengths;

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function RegPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RegPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadReg.
function LoadReg_Callback(hObject, eventdata, handles)
% hObject    handle to LoadReg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname pname] = uigetfile('*.mat','Open a Calibration File');

if fname ~= 0
    CalibFileNames = handles.CalibData.Files;
    if isempty(strfind(pname,handles.CalibData.CalibPath))
        copyfile(fullfile(pname,fname),fullfile(handles.CalibData.CalibPath,fname));
        CalibFileNames{end+1,:} = fname;
        handles.CalibData.Files = CalibFileNames;
        nfiles = length(CalibFileNames);
        set(handles.RegPop,'String',CalibFileNames);
        set(handles.RegPop,'Value',nfiles);
    else
        for i = 1:length(CalibFileNames)
            if strcmp(fname,CalibFileNames{i,:})
                break
            end
        end
        set(handles.RegPop,'Value',i);
    end
    
    
    IN = load(fullfile(pname,fname));
    
    handles.CalibData.RegMat = IN.RegMat;
    handles.CalibData.wavelengths = IN.wavelengths;
    
end
guidata(hObject,handles);

% --- Executes on button press in NewCalibButton.
function NewCalibButton_Callback(hObject, eventdata, handles)
% hObject    handle to NewCalibButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
BeadCalibration_GUI();

curFold = pwd;
cd(handles.CalibData.CalibPath);
CalibFiles = dir('*.mat');

CalibFileNames = cell(size(CalibFiles));

for i = 1:size(CalibFiles)
    CalibFileNames{i,:} = CalibFiles(i).name;
    dates(i) = datenum(CalibFiles(i).date);
end

%Detemine most recent file
mostRecent = find(dates == max(dates));

set(handles.RegPop,'String',CalibFileNames);
set(handles.RegPop,'Value',mostRecent);

%Load in this file
IN = load(CalibFiles(mostRecent).name);

cd(curFold);

handles.CalibData.RegMat = IN.RegMat;
handles.CalibData.wavelengths = IN.wavelengths;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in AlignButton.
function AlignButton_Callback(hObject, eventdata, handles)
% hObject    handle to AlignButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = hObject;

%Determine what wavelengths are present in the calibration file
w_Reg = handles.CalibData.wavelengths;

w_present = unique(w_Reg(:));

for i = 1:length(w_present)
    w_string{i,:} = num2str(w_present(i));
end

%Ask the user to select the wavelengths
w_idx = WavechooseDlg(w_string);

w1 = w_present(w_idx(1));
w2 = w_present(w_idx(2));

%see if the wavelengths are in the same order as the calibration file
RegIndx = find(w_Reg(:,1) == w1 & w_Reg(:,2) == w2,1,'first');
handles.ImagesAlign = handles.Images;
handles.ParticlesAlign = handles.Particles;
handles.TracksAlign = handles.Tracks;
if ~isempty(RegIndx)
    handles.ImagesAlign = handles.Images;
    handles.ParticlesAlign = handles.Particles;
    handles.TracksAlign = handles.Tracks;
    RegMat = handles.CalibData.RegMat{RegIndx,:};
    Roriginal =  imref2d(size(handles.Images.Image1(1).data));
    %Transform the images
    for i = 1:size(handles.Images.Image2,2)
        handles.ImagesAlign.Image2(i).data = imwarp(handles.Images.Image2(i).data,RegMat,'OutputView',Roriginal);
    end
    
    %Transform the particle coordinates
    coord = handles.Particles.Particles2(:,1:2);
    handles.ParticlesAlign.Particles2(:,1:2) = transformPointsForward(RegMat,coord);
    if size(handles.Particles.Particles2,2) > 9
        coord = handles.Particles.Particles2(:,10:11);
        handles.ParticlesAlign.Particles2(:,10:11) = transformPointsForward(RegMat,coord);
    end
    
    %Transform the Track coordinates
    coord = handles.Tracks.Tracks2(:,1:2);
    handles.TracksAlign.Tracks2 = handles.Tracks.Tracks2;
    handles.TracksAlign.Tracks2(:,1:2) = transformPointsForward(RegMat,coord);
else
    RegIndx = find(w_Reg(:,2) == w1 & w_Reg(:,1) == w2,1,'first');
    if isempty(RegIndx)
        errordlg('No mapping found between selected wavelengths');
        return;
    end
    
    handles.ImagesAlign = handles.Images;
    handles.ParticlesAlign = handles.Particles;
    handles.TracksAlign = handles.Tracks;
    RegMat = handles.CalibData.RegMat{RegIndx,:};
    Roriginal =  imref2d(size(handles.Images.Image2(1).data));
    %Transform the images
    for i = 1:size(handles.Images.Image1,2)
        handles.ImagesAlign.Image1(i).data = imwarp(handles.Images.Image1(i).data,RegMat,'OutputView',Roriginal);
    end
    
    %Transform the particle coordinates
        if ~isempty(handles.Particles.Particles1)
        coord = handles.Particles.Particles1(:,1:2);
        handles.ParticlesAlign.Particles1(:,1:2) = transformPointsForward(RegMat,coord);
        if size(handles.Particles.Particles1,2) > 9
            coord = handles.Particles.Particles1(:,10:11);
            handles.ParticlesAlign.Particles1(:,10:11) = transformPointsForward(RegMat,coord);
        end

        %Transform the Track coordinates
        coord = handles.Tracks.Tracks1(:,1:2);
        handles.TracksAlign.Tracks1(:,1:2) = transformPointsForward(RegMat,coord);
        end
    
end

im = merge2colorData(handles.ImagesAlign.Image1(handles.CurrentInd).data,handles.ImagesAlign.Image2(handles.CurrentInd).data,[1,2],handles.clims1,handles.clims2);

axes(handles.TranAxes);
image(im);

% Update handles structure
guidata(hObject, handles);

setappdata(0,'imageStack1',handles.ImagesAlign.Image1);
setappdata(0,'imageStack2',handles.ImagesAlign.Image2);
setappdata(0,'Particles1',handles.ParticlesAlign.Particles1);
setappdata(0,'Particles2',handles.ParticlesAlign.Particles2);
setappdata(0,'Tracks1',handles.TracksAlign.Tracks1);
setappdata(0,'Tracks2',handles.TracksAlign.Tracks2);
    
    


% --- Executes on slider movement.
function ImageSlider_Callback(hObject, eventdata, handles)
% hObject    handle to ImageSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
CurrInd = round(get(hObject,'Value'));

im = merge2colorData(handles.Images.Image1(CurrInd).data,handles.Images.Image2(CurrInd).data,[1 2],handles.clims1,handles.clims2);

axes(handles.OrigAxes);
image(im);

if ~isempty(handles.ImagesAlign)
    im2 = merge2colorData(handles.ImagesAlign.Image1(CurrInd).data,handles.ImagesAlign.Image2(CurrInd).data,[1 2],handles.clims1,handles.clims2);
    axes(handles.TranAxes);
    image(im2);
end
handles.CurrInd = CurrInd;

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ImageSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on button press in DoneButton.
function DoneButton_Callback(hObject, eventdata, handles)
% hObject    handle to DoneButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close();
