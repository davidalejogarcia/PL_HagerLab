function varargout = BeadCalibration_GUI(varargin)
% BEADCALIBRATION_GUI MATLAB code for BeadCalibration_GUI.fig
%      BEADCALIBRATION_GUI, by itself, creates a new BEADCALIBRATION_GUI or raises the existing
%      singleton*.
%
%      H = BEADCALIBRATION_GUI returns the handle to a new BEADCALIBRATION_GUI or the handle to
%      the existing singleton*.
%
%      BEADCALIBRATION_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BEADCALIBRATION_GUI.M with the given input arguments.
%
%      BEADCALIBRATION_GUI('Property','Value',...) creates a new BEADCALIBRATION_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BeadCalibration_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BeadCalibration_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BeadCalibration_GUI

% Last Modified by GUIDE v2.5 08-May-2017 14:16:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BeadCalibration_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @BeadCalibration_GUI_OutputFcn, ...
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


% --- Executes just before BeadCalibration_GUI is made visible.
function BeadCalibration_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BeadCalibration_GUI (see VARARGIN)

% Choose default command line output for BeadCalibration_GUI
handles.output = hObject;
handles.Data.I1 = [];
handles.Data.I2 = [];
handles.Data.I3 = [];
handles.Data.Filepath = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BeadCalibration_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BeadCalibration_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in RunButton.
function RunButton_Callback(hObject, eventdata, handles)
% hObject    handle to RunButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.Data.I1)
    errordlg('No Primary image loaded');
    return;
end
if isempty(handles.Data.I2)
    errordlg('No Secondary image loaded');
    return
end
[movingPts1, fixedPts1] = cpselect(handles.Data.I2,handles.Data.I1,'Wait',true);

%Determine the type of transformation to do
TranSel = get(handles.TranTypePop,'Value');
if TranSel == 1
    TranType = 'affine';
elseif TranSel == 2
    TranType = 'nonreflectivesimilarity';
elseif TranSel == 3
    TranType = 'similarity';
end
RegMat{1,:} = fitgeotrans(movingPts1,fixedPts1,TranType);
% RegMat{1,:} = TwoColorImageRegistration(handles.Data.I1,handles.Data.I2,TranType);
wavelengths(1,:) = [str2double(get(handles.Image1Lambda,'String')), str2double(get(handles.Image2Lambda,'String'))];

% eval(['RegMat_', get(handles.Image1Lambda,'String'), '_', get(handles.Image2Lambda,'String'), ' = TwoColorImageRegistration(handles.Data.I1,handles.Data.I2,TranType);']);

if ~isempty(handles.Data.I3)
    [movingPts2, fixedPts2] = cpselect(handles.Data.I3,handles.Data.I1,movingPts1,fixedPts1,'Wait',true);
    RegMat{2,:} = fitgeotrans(movingPts2,fixedPts2,TranType);
%     RegMat{2,:} = TwoColorImageRegistration(handles.Data.I1,handles.Data.I3,TranType);
%     [movingPts3, fixedPts3] = cpselect(handles.Data.I3,handles.Data.I2,movingPts2,movingPts1,'Wait',true);
    movingPts3 = movingPts2;
    fixedPts3 = movingPts1;
    RegMat{3,:} = fitgeotrans(movingPts3,fixedPts3,TranType);
%     RegMat{3,:} = TwoColorImageRegistration(handles.Data.I2,handles.Data.I3,TranType);
    wavelengths(2,:) = [str2double(get(handles.Image1Lambda,'String')), str2double(get(handles.Image3Lambda,'String'))];
    wavelengths(3,:) = [str2double(get(handles.Image2Lambda,'String')), str2double(get(handles.Image3Lambda,'String'))];

    
    
end
handles.Data.RegMat = RegMat;
handles.Data.wavelengths = wavelengths;
set(handles.SaveButton,'Enable','on');
% Update handles structure
guidata(hObject, handles);



% --- Executes on selection change in TranTypePop.
function TranTypePop_Callback(hObject, eventdata, handles)
% hObject    handle to TranTypePop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns TranTypePop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TranTypePop


% --- Executes during object creation, after setting all properties.
function TranTypePop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TranTypePop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SaveButton.
function SaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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

Im1Lambda = get(handles.Image1Lambda,'String');
Im1Lambda = Im1Lambda;
Im2Lambda = get(handles.Image2Lambda,'String');
Im2Lambda = Im2Lambda;

defName = [CalibPath, filesep, 'Calibration_', date, '_', Im1Lambda,'_',Im2Lambda];
if ~isempty(handles.Data.I3)
    Im3Lambda = get(handles.Image3Lambda,'String');
    Im3Lambda = Im3Lambda;

    defName = [defName, '_', Im3Lambda];
end
defName = [defName,'.mat'];

RegMat = handles.Data.RegMat;

wavelengths = handles.Data.wavelengths;

[fname,pname] = uiputfile('*.mat','Save Calibration Data',defName);

if fname ~= 0
    save([pname, fname],'RegMat','wavelengths');
end


   
    

% --- Executes on button press in Image1Load.
function Image1Load_Callback(hObject, eventdata, handles)
% hObject    handle to Image1Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA
[fname, pname] = uigetfile('*.tif','Open Primary Image, to which other images will be aligned.');

if fname ~= 0
    handles.Data.FilePath = pname;
    handles.Data.I1 = imread([pname, fname]);
    handles.Data.I1 = double(handles.Data.I1);
    handles.Data.I1 = (handles.Data.I1 - min(handles.Data.I1(:)))./(max(handles.Data.I1(:)) - min(handles.Data.I1(:)));
    set(handles.Image1Name,'String',fname);
    set(handles.Image2Load,'Enable','on');
    lambda = inputdlg('Enter the wavelength in nm of this image','Channel Wavelength',1,{'488'});
    set(handles.Image1Lambda,'String',lambda{1});
    set(handles.Image1Lambda,'Enable','on');
    
    
    
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in Image2Load.
function Image2Load_Callback(hObject, eventdata, handles)
% hObject    handle to Image2Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname, pname] = uigetfile('*.tif','Open First Secondary Image',handles.Data.FilePath);

if fname ~= 0
    
    handles.Data.I2 = imread([pname, fname]);
    handles.Data.I2 = double(handles.Data.I2);
    handles.Data.I2 = (handles.Data.I2 - min(handles.Data.I2(:)))./(max(handles.Data.I2(:)) - min(handles.Data.I2(:)));
    set(handles.Image2Name,'String',fname);
    set(handles.Image3Load,'Enable','on');
    lambda = inputdlg('Enter the wavelength in nm of this image','Channel Wavelength',1,{'561'});
    set(handles.Image2Lambda,'String',lambda{1});
    set(handles.Image2Lambda,'Enable','on');
    set(handles.TranTypePop,'Enable','on');
    set(handles.RunButton,'Enable','on');
    
    
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in Image3Load.
function Image3Load_Callback(hObject, eventdata, handles)
% hObject    handle to Image3Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname, pname] = uigetfile('*.tif','Open Other Secondary Image',handles.Data.FilePath);

if fname ~= 0
    
    handles.Data.I3 = imread([pname, fname]);
    handles.Data.I3 = double(handles.Data.I3);
    handles.Data.I3 = (handles.Data.I3 - min(handles.Data.I3(:)))./(max(handles.Data.I3(:)) - min(handles.Data.I3(:)));
    set(handles.Image3Name,'String',fname);
%     set(handles.Image3Load,'Enable','on');
    lambda = inputdlg('Enter the wavelength in nm of this image','Channel Wavelength',1,{'647'});
    set(handles.Image3Lambda,'String',lambda{1});
    set(handles.Image3Lambda,'Enable','on');
    set(handles.TranTypePop,'Enable','on');
    set(handles.RunButton,'Enable','on');
    
    
end

% Update handles structure
guidata(hObject, handles);


function Image1Lambda_Callback(hObject, eventdata, handles)
% hObject    handle to Image1Lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Image1Lambda as text
%        str2double(get(hObject,'String')) returns contents of Image1Lambda as a double


% --- Executes during object creation, after setting all properties.
function Image1Lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Image1Lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Image2Lambda_Callback(hObject, eventdata, handles)
% hObject    handle to Image2Lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Image2Lambda as text
%        str2double(get(hObject,'String')) returns contents of Image2Lambda as a double


% --- Executes during object creation, after setting all properties.
function Image2Lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Image2Lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Image3Lambda_Callback(hObject, eventdata, handles)
% hObject    handle to Image3Lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Image3Lambda as text
%        str2double(get(hObject,'String')) returns contents of Image3Lambda as a double


% --- Executes during object creation, after setting all properties.
function Image3Lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Image3Lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadRaw.
function LoadRaw_Callback(hObject, eventdata, handles)
% hObject    handle to LoadRaw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname, pname] = uigetfile('*.tif','Open Image Stack');
if fname ~= 0
    %re-initialize buttons and displays
    set(handles.RunButton,'Enable','off');
    set(handles.Image2Load,'Enable','off');
    set(handles.Image3Load,'Enable','off');
    set(handles.SaveButton,'Enable','off');
    set(handles.TranTypePop,'Enable','off');
    
    set(handles.Image1Name,'String','');
    set(handles.Image2Name,'String','');
    set(handles.Image3Name,'String','');
    
    set(handles.Image1Lambda,'String','','Enable','off');
    set(handles.Image2Lambda,'String','','Enable','off');
    set(handles.Image3Lambda,'String','','Enable','off');
    drawnow;
    
    %save the image name and path to handles structure
    handles.Data.ImageFilename = fname;
    handles.Data.ImageFilepath = pname;
    
    
    channelOrderStrCell = inputdlg('Enter the channel order for the loaded image stack as comma seperataed values','Channel Order',1,{'647,561,488'});
    [ImageStack, nImages] = TIFread(fullfile(pname,fname));
    handles.Data.RawImages = ImageStack;
    channelOrderStr = channelOrderStrCell{1,:};
    %Determine how many channels are specified
    commaPos = strfind(channelOrderStr,',');
    fname2 = cell(3,1);
    if length(commaPos) == 1
        nChan = 2;
        wl1_num = str2double(channelOrderStr(1:commaPos(1)-1));
        wl1_str = num2str(wl1_num);
        fname2{1,:} = [handles.Data.ImageFilename(1:end-4),'_',wl1_str,'.tif'];
        wl2_num = str2double(channelOrderStr(commaPos(1)+1:end));
        wl2_str = num2str(wl2_num);
        fname2{2,:} = [handles.Data.ImageFilename(1:end-4),'_',wl2_str,'.tif'];
        
        handles.Data.wls = [wl1_num;wl2_num];
        
    elseif length(commaPos) == 2
        nChan = 3;
        wl1_num = str2double(channelOrderStr(1:commaPos(1)-1));
        wl1_str = num2str(wl1_num);
        fname2{1,:} = [handles.Data.ImageFilename(1:end-4),'_',wl1_str,'.tif'];
        wl2_num = str2double(channelOrderStr(commaPos(1)+1:commaPos(2)-1));
        wl2_str = num2str(wl2_num);
        fname2{2,:} = [handles.Data.ImageFilename(1:end-4),'_',wl2_str,'.tif'];
        wl3_num = str2double(channelOrderStr(commaPos(2)+1:end));
        wl3_str = num2str(wl3_num);
        fname2{3,:} = [handles.Data.ImageFilename(1:end-4),'_',wl3_str,'.tif'];
        
        handles.Data.wls = [wl1_num;wl2_num;wl3_num];
    elseif isempty(commaPos)
        errordlg('A single channel was specified, or commas were missing','Unrecognized string format');
    elseif length(commaPos) > 2
        errordlg('Splitting can currently be done for a maximum of 3 wavelengths ','Too many input channels');
    end
    
    %Initialize image files
    handles.Data.ImAll = cell(nChan,1);
    counter = ones(nChan,1);
    matlabpath = path;
    sc_pos = strfind(matlabpath,pathsep);
    
    handles.Software.TempSaveDir = [matlabpath(1:sc_pos(1)-1), filesep, 'tmp',filesep]; %temporary save directory to avoid the error when saving to a folder that is open in Windows explorer
    TempDirExists = dir(handles.Software.TempSaveDir);
    if isempty(TempDirExists)
        mkdir(handles.Software.TempSaveDir);
    end
    
    for i = 1:nImages
        im_tmp = handles.Data.RawImages(i).data;
        if rem(i,nChan) == 1
            
            handles.Data.I1(:,:,counter(1)) = im_tmp;
            counter(1) = counter(1) + 1;
            fnum = 1;
        elseif rem(i,nChan) == 2
            
            handles.Data.I2(:,:,counter(2)) = im_tmp;
            counter(2) = counter(2) + 1;
            fnum = 2;
        elseif rem(i,nChan) == 0 && nChan == 2
            
            handles.Data.I2(:,:,counter(2)) = im_tmp;
            counter(2) = counter(2) + 1;
            fnum = 2;
        elseif rem(i,nChan) == 0 && nChan == 3
            
            handles.Data.I3(:,:,counter(3)) = im_tmp;
            counter(3) = counter(3) + 1;
            fnum = 3;
        end
        %     set(handles.StatusText,'String',['Writing Channel data: ',num2str(handles.Data.wls(fnum)),'...']);
        %     drawnow
        %     imwrite(im_tmp,[handles.Software.TempSaveDir, fname2{fnum,:}],'WriteMode','append');
        %     set(handles.StatusText,'String',['Writing Channel data: ',num2str(handles.Data.wls(fnum)),'...Done']);
        %     drawnow
    end
    % Create the images to be used for registration (either the only image
    % or a max. projection).
    if nImages == nChan
        handles.Data.I1 = handles.Data.I1(:,:,1);
        
        handles.Data.I2 = handles.Data.I2(:,:,1);
        
        if nChan > 2
            handles.Data.I3 = handles.Data.I3(:,:,1);
            
        end
    else
        handles.Data.I1 = max(handles.Data.I1,[],3);
        handles.Data.I2 = max(handles.Data.I2,[],3);
        if nChan > 2
            handles.Data.I3 = max(handles.Data.I3,[],3);
        end
    end
    handles.Data.I1 = double(handles.Data.I1);
    handles.Data.I1 = (handles.Data.I1 - min(handles.Data.I1(:)))./(max(handles.Data.I1(:)) - min(handles.Data.I1(:)));
    handles.Data.I2 = double(handles.Data.I2);
    handles.Data.I2 = (handles.Data.I2 - min(handles.Data.I2(:)))./(max(handles.Data.I2(:)) - min(handles.Data.I2(:)));
    if nChan > 2
        handles.Data.I3 = double(handles.Data.I3);
        handles.Data.I3 = (handles.Data.I3 - min(handles.Data.I3(:)))./(max(handles.Data.I3(:)) - min(handles.Data.I3(:)));
    end
    
    set(handles.Image1Lambda,'String',num2str(handles.Data.wls(1)));
    set(handles.Image1Name,'String',[fname, ' Chan 1'],'Enable','on');
    set(handles.Image1Load,'Enable','on');
    set(handles.Image2Lambda,'String',num2str(handles.Data.wls(2)));
    set(handles.Image2Name,'String',[fname, ' Chan 2'],'Enable','on');
    set(handles.Image2Load,'Enable','on');
    set(handles.Image3Load,'Enable','on');
    if nChan > 2
        set(handles.Image3Lambda,'String',num2str(handles.Data.wls(3)));
        set(handles.Image3Name,'String',[fname, ' Chan 3'],'Enable','on');
    end
    set(handles.RunButton,'Enable','on');
    set(handles.TranTypePop,'Enable','on');
    guidata(hObject,handles);
    
end
