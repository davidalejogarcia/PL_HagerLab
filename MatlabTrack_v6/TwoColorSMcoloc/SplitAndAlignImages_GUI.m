function varargout = SplitAndAlignImages_GUI(varargin)
% SPLITANDALIGNIMAGES_GUI MATLAB code for SplitAndAlignImages_GUI.fig
%      SPLITANDALIGNIMAGES_GUI, by itself, creates a new SPLITANDALIGNIMAGES_GUI or raises the existing
%      singleton*.
%
%      H = SPLITANDALIGNIMAGES_GUI returns the handle to a new SPLITANDALIGNIMAGES_GUI or the handle to
%      the existing singleton*.
%
%      SPLITANDALIGNIMAGES_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPLITANDALIGNIMAGES_GUI.M with the given input arguments.
%
%      SPLITANDALIGNIMAGES_GUI('Property','Value',...) creates a new SPLITANDALIGNIMAGES_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SplitAndAlignImages_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SplitAndAlignImages_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SplitAndAlignImages_GUI

% Last Modified by GUIDE v2.5 05-May-2017 14:36:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SplitAndAlignImages_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @SplitAndAlignImages_GUI_OutputFcn, ...
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


% --- Executes just before SplitAndAlignImages_GUI is made visible.
function SplitAndAlignImages_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SplitAndAlignImages_GUI (see VARARGIN)

% Choose default command line output for SplitAndAlignImages_GUI
handles.Reg.file = '';
handles.Data.ImageFilename = '';
handles.Data.ImAlign = [];
handles.Data.ImAll = [];
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SplitAndAlignImages_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SplitAndAlignImages_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ImageFileSelectButton.
function ImageFileSelectButton_Callback(hObject, eventdata, handles)
% hObject    handle to ImageFileSelectButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fname, pname] = uigetfile('*.tif','Open Image Stack');

if fname ~= 0
    %save the image name and path to handles structure
    handles.Data.ImageFilename = fname;
    handles.Data.ImageFilepath = pname;
    %read in the images
    statusStr = ['Loading image: ',fullfile(pname,fname),'...'];
    set(handles.StatusText,'String',statusStr);
    drawnow
    [ImageStack, nImages] = TIFread(fullfile(pname,fname));
    set(handles.StatusText,'String',['Loading image: ',fullfile(pname,fname),'...Done']);
    drawnow
    handles.Data.RawImages = ImageStack;
    handles.Data.nImages = nImages;
    
    %add the file name to the text box & enable the split colors button
    set(handles.ImageFileText,'String',fullfile(pname,fname));
    set(handles.SplitColorsButton,'Enable','on');
    
    if ~strcmp(handles.Reg.file,'')
        set(handles.SplitRegisterButton,'Enable','on');
    end
    
    guidata(hObject,handles);
end


function ImageFileText_Callback(hObject, eventdata, handles)
% hObject    handle to ImageFileText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ImageFileText as text
%        str2double(get(hObject,'String')) returns contents of ImageFileText as a double

if ~isempty(get(hObject,'String'))
    set(handles.StatusText,'String',['Loading image: ',get(hObject,'String'),'...']);
    drawnow
    
    try [ImageStack, nImages] = TIFread(get(hObject,'String'));
        
    catch
        errordlg('File does not exist or is not an image','Error opening file');
    end
    set(handles.StatusText,'String',['Loading image: ',get(hObject,'String'),'...Done']);
    drawnow
    handles.Data.RawImages = ImageStack;
    handles.Data.nImages = nImages;
    %add the file name to the text box & enable the split colors button
    
    set(handles.SplitColorsButton,'Enable','on');
    
    %Extract file name & path
    file_full = get(hObject,'String');
    slash_pos = strfind(file_full,filesep);
    
    fname = file_full(slash_pos(end+1):end);
    pname = file_full(1:slash_pos(end-1));
    %save the image name and path to handles structure
    handles.Data.ImageFilename = fname;
    handles.Data.ImageFilepath = pname;
    
    if ~strcmp(handles.Reg.file,'')
        set(handles.SplitRegisterButton,'Enable','on');
    end
    
    guidata(hObject,handles);
    
end

% --- Executes during object creation, after setting all properties.
function ImageFileText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageFileText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CalibFileSelectButton.
function CalibFileSelectButton_Callback(hObject, eventdata, handles)
% hObject    handle to CalibFileSelectButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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
[fname, pname] = uigetfile(fullfile(CalibPath,'*.mat'),'Open Mat file with Registration info');

if fname ~= 0
    handles.Reg.file = fullfile(pname,fname);
    set(handles.CalibFileText,'String',handles.Reg.file);
    set(handles.StatusText,'String',['Loading Alignment data: ',fullfile(pname,fname),'...']);
    drawnow
    RegData = load(fullfile(pname,fname));
    set(handles.StatusText,'String',['Loading Alignment data: ',fullfile(pname,fname),'...Done']);
    drawnow
    
    handles.Reg.RegMat = RegData.RegMat;
    handles.Reg.wls = RegData.wavelengths;
    
    %     if ~isempty(handles.Data.ImAlign)
    set(handles.RegisterButton,'Enable','on');
    %     end
    
    if ~strcmp(handles.Data.ImageFilename,'')
        set(handles.SplitRegisterButton,'Enable','on');
    end
    guidata(hObject,handles);
end


function CalibFileText_Callback(hObject, eventdata, handles)
% hObject    handle to CalibFileText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CalibFileText as text
%        str2double(get(hObject,'String')) returns contents of CalibFileText as a double


% --- Executes during object creation, after setting all properties.
function CalibFileText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CalibFileText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SplitColorsButton.
function SplitColorsButton_Callback(hObject, eventdata, handles)
% hObject    handle to SplitColorsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
channelOrderStrCell = inputdlg('Enter the channel order for the loaded image stack as comma seperataed values','Channel Order',1,{'647,561,488'});

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

for i = 1:handles.Data.nImages
    im_tmp = handles.Data.RawImages(i).data;
    if rem(i,nChan) == 1
        
        handles.Data.ImAll{1}(:,:,counter(1)) = im_tmp;
        counter(1) = counter(1) + 1;
        fnum = 1;
    elseif rem(i,nChan) == 2
        
        handles.Data.ImAll{2}(:,:,counter(2)) = im_tmp;
        counter(2) = counter(2) + 1;
        fnum = 2;
    elseif rem(i,nChan) == 0 && nChan == 2
        
        handles.Data.ImAll{2}(:,:,counter(2)) = im_tmp;
        counter(2) = counter(2) + 1;
        fnum = 2;
    elseif rem(i,nChan) == 0 && nChan == 3
        
        handles.Data.ImAll{3}(:,:,counter(3)) = im_tmp;
        counter(3) = counter(3) + 1;
        fnum = 3;
    end
    %     set(handles.StatusText,'String',['Writing Channel data: ',num2str(handles.Data.wls(fnum)),'...']);
    %     drawnow
    %     imwrite(im_tmp,[handles.Software.TempSaveDir, fname2{fnum,:}],'WriteMode','append');
    %     set(handles.StatusText,'String',['Writing Channel data: ',num2str(handles.Data.wls(fnum)),'...Done']);
    %     drawnow
end
if get(handles.CutFirstFrame,'Value') == 1
    frame1 = 2;
else
    frame1 = 1;
end

for i = 1:nChan
    set(handles.StatusText,'String',['Writing Channel ', num2str(i),' data...']);
    drawnow
    for j = frame1:counter(1)-1
        im_tmp2 = handles.Data.ImAll{i}(:,:,j);
        imwrite(im_tmp2,[handles.Software.TempSaveDir, fname2{i,:}],'WriteMode','append');
    end
    
    movefile([handles.Software.TempSaveDir, filesep, fname2{i,:}],[handles.Data.ImageFilepath, filesep, fname2{i,:}]);
    set(handles.StatusText,'String',['Writing Channel ', num2str(i), ' data...Done']);
    drawnow
end
set(handles.StatusText,'String','Channel splitting Done!');
guidata(hObject,handles);

% --- Executes on button press in RegisterButton.
function RegisterButton_Callback(hObject, eventdata, handles)
% hObject    handle to RegisterButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Check if there is some image data from the current session
if ~isempty(handles.Data.ImAll)
    CurOrOldChoice = questdlg('Do you want to perform alignment on the images just split?','What do you want to align?','Current','Select','Current');
else
    CurOrOldChoice = 'Select';
end

if ~isfield(handles.Data,'ImageFilepath') || isempty(handles.Data.ImageFilepath)
    handles.Data.ImageFilepath = pwd;
    handles.Data.ImageFilename = [];
end
if size(handles.Reg.wls,1) == 3
    wl = [handles.Reg.wls(1,1);handles.Reg.wls(1,2);handles.Reg.wls(2,2)];
elseif size(handles.Reg.wls,1) == 1
    wl = [handles.Reg.wls(1,1);handles.Reg.wls(1,2)];
end
pname = '';
wl_indx = 1;
if strcmp(CurOrOldChoice,'Select')
    for i = 1:length(wl)
        openFileStr = ['Select file for ',num2str(wl(i)), ' channel'];
        if strcmp(pname,'')
            [fname, pname] = uigetfile('*.tif',openFileStr,handles.Data.ImageFilepath);
        else
            [fname, pname] = uigetfile('*.tif',openFileStr,pname);
        end
        if fname ~= 0
            handles.Data.ImageFilepath = pname;
            handles.Data.ImageFilename = fname;
            wluse(wl_indx,1) = wl(i);
            [imStack,nImages] = TIFread([pname,filesep,fname]);
            for j = 1:nImages
                Im2Align{wl_indx,:}(:,:,j) = imStack(j).data;
            end
            continueFlag = 1;
            wl_indx = wl_indx + 1;
        elseif i == length(wl) && isempty(handles.Data.ImageFilename)
            continueFlag = 0;
        else
            continueFlag = 1;
            pname = handles.Data.ImageFilepath;
        end
    end
else
    pname = handles.Data.ImageFilepath;
    Im2Align = handles.Data.ImAll;
    wluse = handles.Data.wls;
    continueFlag = 1;
end


%align the images
if continueFlag == 1
    matlabpath = path;
    sc_pos = strfind(matlabpath,pathsep);

    handles.Software.TempSaveDir = [matlabpath(1:sc_pos(1)-1), filesep, 'tmp',filesep]; %temporary save directory to avoid the error when saving to a folder that is open in Windows explorer
    TempDirExists = dir(handles.Software.TempSaveDir);
    if isempty(TempDirExists)
        mkdir(handles.Software.TempSaveDir);
    end
%     mkdir([pname, filesep,'aligned']);
    for i = 1:size(wluse,1)
        set(handles.StatusText,'String',['Aligning Channel ', num2str(i), 'data...']);
        drawnow
        if wluse(i) ~= wl(1);
            Reg_ind = find(handles.Reg.wls(:,1) == wl(1) & handles.Reg.wls(:,2) == wluse(i));
            RegMat = handles.Reg.RegMat{Reg_ind,:};
            Roriginal = imref2d([size(Im2Align{1},1),size(Im2Align{1},2)]);
            for j = 1:size(Im2Align{i},3)
                handles.Data.ImAlign{i}(:,:,j) = imwarp(Im2Align{i}(:,:,j),RegMat,'OutputView',Roriginal);
                imwrite(handles.Data.ImAlign{i}(:,:,j),[handles.Software.TempSaveDir,filesep,handles.Data.ImageFilename(1:end-4),'_',num2str(wluse(i)),'_aligned.tif'],'WriteMode','append');
            end
            movefile([handles.Software.TempSaveDir, handles.Data.ImageFilename(1:end-4),'_',num2str(wluse(i)),'_aligned.tif'],...
                [handles.Data.ImageFilepath, filesep, handles.Data.ImageFilename(1:end-4),'_',num2str(wluse(i)),'_aligned.tif']);
%         else
%             copyfile([pname,filesep,handles.Data.ImageFilename(1:end-4),'_',num2str(wluse(i)),'.tif'],...
%                 [pname,filesep,'aligned',filesep,handles.Data.ImageFilename(1:end-4),'_',num2str(wluse(i)),'.tif']);
        end
        set(handles.StatusText,'String',['Aligning Channel ', num2str(i), 'data...Done']);
        drawnow;
    end
end
set(handles.StatusText,'String','Channel Alignment Done!');
guidata(hObject,handles);

% --- Executes on button press in SplitRegisterButton.
function SplitRegisterButton_Callback(hObject, eventdata, handles)
% hObject    handle to SplitRegisterButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
channelOrderStrCell = inputdlg('Enter the channel order for the loaded image stack as comma seperataed values','Channel Order',1,{'647,561,488'});
wluse = 0;
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

wluse = handles.Data.wls;

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
set(handles.StatusText,'String','Splitting Channels...');
drawnow
for i = 1:handles.Data.nImages
    im_tmp = handles.Data.RawImages(i).data;
    if rem(i,nChan) == 1
        
        handles.Data.ImAll{1}(:,:,counter(1)) = im_tmp;
        counter(1) = counter(1) + 1;
        fnum = 1;
    elseif rem(i,nChan) == 2
        
        handles.Data.ImAll{2}(:,:,counter(2)) = im_tmp;
        counter(2) = counter(2) + 1;
        fnum = 2;
    elseif rem(i,nChan) == 0 && nChan == 2
        
        handles.Data.ImAll{2}(:,:,counter(2)) = im_tmp;
        counter(2) = counter(2) + 1;
        fnum = 2;
    elseif rem(i,nChan) == 0 && nChan == 3
        
        handles.Data.ImAll{3}(:,:,counter(3)) = im_tmp;
        counter(3) = counter(3) + 1;
        fnum = 3;
    end
    
%     imwrite(im_tmp,[handles.Software.TempSaveDir, fname2{fnum}],'WriteMode','append');
end
Im2Align = handles.Data.ImAll;
set(handles.StatusText,'String','Splitting Channels...Done');
drawnow

if get(handles.CutFirstFrame,'Value') == 1
    frame1 = 2;
else
    frame1 = 1;
end
for i = 1:nChan
    set(handles.StatusText,'String',['Writing Channel ', num2str(i), ' data...']);
    drawnow
    for j = frame1:counter(1)-1
        im_tmp2 = handles.Data.ImAll{i}(:,:,j);
        imwrite(im_tmp2,[handles.Software.TempSaveDir, fname2{i,:}],'WriteMode','append');
    end
    
    movefile([handles.Software.TempSaveDir, filesep, fname2{i,:}],[handles.Data.ImageFilepath, filesep, fname2{i,:}]);
    im_max = max(handles.Data.ImAll{i},[],3);
    imwrite(im_max,[handles.Data.ImageFilepath, filesep, 'MAX_',fname2{i,:}]);
%     im_sum = sum(handles.Data.ImAll{i},3);
%     imwrite((im_sum),[handles.Data.ImageFilepath, filesep, 'SUM_',fname2{i,:}]);
    im_avg = mean(handles.Data.ImAll{i},3);
    imwrite(uint16(im_avg),[handles.Data.ImageFilepath, filesep, 'AVG_',fname2{i,:}]);
    
    set(handles.StatusText,'String',['Writing Channel ', num2str(i), ' data...Done']);
    drawnow
end

% mkdir([handles.Data.ImageFilepath, filesep,'aligned']);
if size(handles.Reg.wls,1) == 3
    wl = [handles.Reg.wls(1,1);handles.Reg.wls(1,2);handles.Reg.wls(2,2)];
elseif size(handles.Reg.wls,1) == 1
    wl = [handles.Reg.wls(1,1);handles.Reg.wls(1,2)];
end

for i = 1:size(wluse,1)
    set(handles.StatusText,'String',['Aligning Channel ', num2str(i), ' data...']);
    drawnow
    if wluse(i) ~= wl(1);
        
        Reg_ind = find(handles.Reg.wls(:,1) == wl(1) & handles.Reg.wls(:,2) == wluse(i));
        RegMat = handles.Reg.RegMat{Reg_ind,:};
        Roriginal = imref2d([size(Im2Align{1},1),size(Im2Align{1},2)]);
        for j = frame1:size(Im2Align{i},3)
            handles.Data.ImAlign{i}(:,:,j) = imwarp(Im2Align{i}(:,:,j),RegMat,'OutputView',Roriginal);
            imwrite(handles.Data.ImAlign{i}(:,:,j),[handles.Software.TempSaveDir,filesep,handles.Data.ImageFilename(1:end-4),'_',num2str(wluse(i)),'_aligned.tif'],'WriteMode','append');
        end
        movefile([handles.Software.TempSaveDir, handles.Data.ImageFilename(1:end-4),'_',num2str(wluse(i)),'_aligned.tif'],[handles.Data.ImageFilepath, filesep, handles.Data.ImageFilename(1:end-4),'_',num2str(wluse(i)),'_aligned.tif']);
        im_max = max(handles.Data.ImAlign{i},[],3);
        imwrite(im_max,[handles.Data.ImageFilepath, filesep, 'MAX_',handles.Data.ImageFilename(1:end-4),'_',num2str(wluse(i)),'_aligned.tif']);
%         im_sum = sum(handles.Data.ImAlign{i},3);
%         imwrite(im_sum,[handles.Data.ImageFilepath, filesep, 'SUM_',handles.Data.ImageFilename(1:end-4),'_',num2str(wluse(i)),'_aligned.tif']);
        im_avg = mean(handles.Data.ImAlign{i},3);
        imwrite(uint16(im_avg),[handles.Data.ImageFilepath, filesep, 'AVG_',handles.Data.ImageFilename(1:end-4),'_',num2str(wluse(i)),'_aligned.tif']);
        set(handles.StatusText,'String',['Aligning Channel ', num2str(i), 'data...Done']);
%     else
%         copyfile([handles.Data.ImageFilepath,handles.Data.ImageFilename(1:end-4),'_',num2str(wluse(i)),'.tif'],...
%             [handles.Data.ImageFilepath,filesep,filesep,handles.Data.ImageFilename(1:end-4),'_',num2str(wluse(i)),'_aligned.tif'])
    end
     set(handles.StatusText,'String',['Aligning Channel ', num2str(i), ' data...Done']);
end
set(handles.StatusText,'String','Splitting and Aligning Done!');
guidata(hObject,handles);
% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function NewCalibMenu_Callback(hObject, eventdata, handles)
% hObject    handle to NewCalibMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

BeadCalibration_GUI;


% --------------------------------------------------------------------
function LoadImageMenu_Callback(hObject, eventdata, handles)
% hObject    handle to LoadImageMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname, pname] = uigetfile('*.tif','Open Image Stack');

if fname ~= 0
    %save the image name and path to handles structure
    handles.Data.ImageFilename = fname;
    handles.Data.ImageFilepath = pname;
    %read in the images
    set(handles.StatusText,'String',['Loading image: ',fullfile(pname,fname),'...']);
    drawnow
    [ImageStack, nImages] = TIFread(fullfile(pname,fname));
    set(handles.StatusText,'String',['Loading image: ',fullfile(pname,fname),'...Done']);
    drawnow
    handles.Data.RawImages = ImageStack;
    handles.Data.nImages = nImages;
    
    %add the file name to the text box & enable the split colors button
    set(handles.ImageFileText,'String',fullfile(pname,fname));
    set(handles.SplitColorsButton,'Enable','on');
    
    if ~strcmp(handles.Reg.file,'')
        set(handles.SplitRegisterButton,'Enable','on');
    end
    
    guidata(hObject,handles);
end

% --------------------------------------------------------------------
function LoadCalibMenu_Callback(hObject, eventdata, handles)
% hObject    handle to LoadCalibMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname, pname] = uigetfile('*.mat','Open Mat file with Registration info');

if fname ~= 0
    handles.Reg.file = fullfile(pname,fname);
    set(handles.StatusText,'String',['Loading Alignment data: ',fullfile(pname,fname),'...']);
    drawnow
    RegData = load(fullfile(pname,fname));
    set(handles.StatusText,'String',['Loading Alignment data: ',fullfile(pname,fname),'...Done']);
    drawnow
    handles.Reg.RegMat = RegData.RegMat;
    handles.Reg.wls = RegData.wavelengths;
    
    if ~isempty(handles.Data.ImAlign)
        set(handles.RegisterButton,'Enable','on');
    end
    
    if ~strcmp(handles.Data.ImageFileName,'')
        set(handles.SplitRegisterButton,'Enable','on');
    end
    
end


% --------------------------------------------------------------------
function BatchMenu_Callback(hObject, eventdata, handles)
% hObject    handle to BatchMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%Ask if we should use the currently loaded Calibration file
if ~strcmp(handles.Reg.file,'')
    CurRegChoose = questdlg(['Do you want to use the loaded Calibration file: ',handles.Reg.file,' or a different one?'],'Calibration file','Current','Select','Current');
else
    CurRegChoose = 'Select';
end

CutFirstChoose = questdlg('Do you want to remove the first frame from each movie?', 'Cut First Frame','Yes','No','No');

if strcmp(CutFirstChoose,'Yes')
    frame1 = 2;
else
    frame1 = 1;
end

if strcmp(CurRegChoose,'Select') %if the user wants to use a different file, open the dialog selection tool and load in
    [fname,pname] = uigetfile('*.mat','Open Calibration File');
    if fname ~= 0
        set(handles.StatusText,'String',['Loading Alignment data: ',fullfile(pname,fname),'...']);
        drawnow
        RegData = load(fullfile(pname,fname));
        RegMat = RegData.RegMat;
        Reg_wls = RegData.wavelengths;
        set(handles.StatusText,'String',['Loading Alignment data: ',fullfile(pname,fname),'...Done']);
        drawnow
        continueFlag = 1;
    else
        continueFlag = 0;
        set(handles.StatusText,'String','Batch Processing Aborted');
        drawnow
    end
    
else
    RegMat = handles.Reg.RegMat;
    Reg_wls = handles.Reg.wls;
    continueFlag = 1;
end

if continueFlag == 1
    if size(Reg_wls,1) == 3
        wl = [Reg_wls(1,1);Reg_wls(1,2);Reg_wls(2,2)];
    elseif size(Reg_wls,1) == 1
        wl = [Reg_wls(1,1);Reg_wls(1,2)];
    end
%     wl = [Reg_wls(1,1);Reg_wls(1,2);Reg_wls(2,2)];
    FileRootPath = uigetdir(pwd,'Select Root Folder containing images');
    if FileRootPath ~= 0
        set(handles.StatusText,'String','Recusively searching Root Folder for images...');
        drawnow
        Dcur = pwd;
        Folds2Search = {FileRootPath};
        InternalFolds = {FileRootPath};
        %Recursively search for possible images to load
        while ~isempty(InternalFolds)
            InternalFolds_new = {};
            fold_ind = 1;
            for i = 1:size(InternalFolds)
                cd(InternalFolds{i,:});
                
                FoldCont = dir;
                FoldCont(~[FoldCont.isdir]) = [];
                tf = ismember({FoldCont.name},{'.','..'});
                FoldCont(tf) = [];
                if ~isempty(FoldCont)
                    for j = 1:size(FoldCont,1)
                        InternalFolds_new = [InternalFolds_new; [InternalFolds{i,:}, filesep, FoldCont(j).name]];
                    end
                end
            end
            InternalFolds = InternalFolds_new;
            Folds2Search = [Folds2Search; InternalFolds];
        end
        IMpaths = {};
        IMfiles = {};
        for i = 1:size(Folds2Search,1)
            
            cd(Folds2Search{i,:})
            imsPresent = dir('*.tif');
            if ~isempty(imsPresent)
                for j = 1:size(imsPresent,1)
                    IMpaths = [IMpaths; Folds2Search(i)];
                    IMfiles = [IMfiles; imsPresent(j).name];
                end
            end
        end
        
        
        set(handles.StatusText,'String',['Found ', num2str(size(IMfiles,1)), ' images']);
        drawnow
        
        channelOrderStrCell = inputdlg('Enter the channel order for the loaded image stack as comma seperataed values','Channel Order',1,{'647,561,488'});
        
        channelOrderStr = channelOrderStrCell{1,:};
        %Determine how many channels are specified
        commaPos = strfind(channelOrderStr,',');
        
        if length(commaPos) == 1
            nChan = 2;
            wl1_num = str2double(channelOrderStr(1:commaPos(1)-1));
            wl1_str = num2str(wl1_num);
            
            wl2_num = str2double(channelOrderStr(commaPos(1)+1:end));
            wl2_str = num2str(wl2_num);
            
            
            wluse = [wl1_num;wl2_num];
            
        elseif length(commaPos) == 2
            nChan = 3;
            wl1_num = str2double(channelOrderStr(1:commaPos(1)-1));
            wl1_str = num2str(wl1_num);
           
            wl2_num = str2double(channelOrderStr(commaPos(1)+1:commaPos(2)-1));
            wl2_str = num2str(wl2_num);
            
            wl3_num = str2double(channelOrderStr(commaPos(2)+1:end));
            wl3_str = num2str(wl3_num);
            
            
            wluse = [wl1_num;wl2_num;wl3_num];
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
        for ImIx = 1:size(IMfiles,1)
            clear ImAll;
            clear ImAlign;
            counter = ones(nChan,1);
            set(handles.StatusText,'String',['Reading Image ', num2str(ImIx), ' of ', num2str(size(IMfiles,1)),'...']);
            drawnow
            [imageStack,nImages] = TIFread([IMpaths{ImIx,:},filesep,IMfiles{ImIx,:}]);
            set(handles.StatusText,'String',['Reading Image ', num2str(ImIx), ' of ', num2str(size(IMfiles,1)),'...Done']);
            drawnow
            fname2 = cell(3,1);
            if length(commaPos) == 1
                
                fname2{1,:} = [IMfiles{ImIx,:}(1:end-4),'_',wl1_str,'.tif'];
                
                fname2{2,:} = [IMfiles{ImIx,:}(1:end-4),'_',wl2_str,'.tif'];
                
                
                
            elseif length(commaPos) == 2
                
                fname2{1,:} = [IMfiles{ImIx,:}(1:end-4),'_',wl1_str,'.tif'];
                
                fname2{2,:} = [IMfiles{ImIx,:}(1:end-4),'_',wl2_str,'.tif'];
                
                fname2{3,:} = [IMfiles{ImIx,:}(1:end-4),'_',wl3_str,'.tif'];
                
                
            end
            set(handles.StatusText,'String',['Splitting Image ', num2str(ImIx), ' of ', num2str(size(IMfiles,1)),'...']);
            drawnow
            for i = 1:nImages
                im_tmp = imageStack(i).data;
                if rem(i,nChan) == 1
                    
                    ImAll{1}(:,:,counter(1)) = im_tmp;
                    counter(1) = counter(1) + 1;
                    fnum = 1;
                elseif rem(i,nChan) == 2
                    
                    ImAll{2}(:,:,counter(2)) = im_tmp;
                    counter(2) = counter(2) + 1;
                    fnum = 2;
                elseif rem(i,nChan) == 0 && nChan == 2
                    
                    ImAll{2}(:,:,counter(2)) = im_tmp;
                    counter(2) = counter(2) + 1;
                    fnum = 2;
                elseif rem(i,nChan) == 0 && nChan == 3
                    
                    ImAll{3}(:,:,counter(3)) = im_tmp;
                    counter(3) = counter(3) + 1;
                    fnum = 3;
                end
                
                
            end
            
            for i = 1:nChan
                for j = frame1:counter(i)-1
                    im_tmp2 = ImAll{i}(:,:,j);
                    imwrite(im_tmp2,[handles.Software.TempSaveDir, fname2{i,:}],'WriteMode','append');
                end
                
                movefile([handles.Software.TempSaveDir, filesep, fname2{i,:}],[IMpaths{ImIx,:}, filesep, fname2{i,:}]);
            end
            Im2Align = ImAll;
            set(handles.StatusText,'String',['Aligning Image ', num2str(ImIx), ' of ', num2str(size(IMfiles,1)),'...']);
            drawnow
%             mkdir([IMpaths{ImIx,:}, filesep,'aligned']);
            for i = 1:size(wluse,1)
                if wluse(i) ~= Reg_wls(1);
                    Reg_ind = find(Reg_wls(:,1) == wl(1) & Reg_wls(:,2) == wluse(i));
                    RegMat_use = RegMat{Reg_ind,:};
                    Roriginal = imref2d([size(Im2Align{1},1),size(Im2Align{1},2)]);
                    for j = frame1:size(Im2Align{i},3)
                        ImAlign{i}(:,:,j) = imwarp(Im2Align{i}(:,:,j),RegMat_use,'OutputView',Roriginal);
                        imwrite(ImAlign{i}(:,:,j),[handles.Software.TempSaveDir,IMfiles{ImIx,:}(1:end-4),'_',num2str(wluse(i)),'_aligned.tif'],'WriteMode','append');
                    end
                    movefile([handles.Software.TempSaveDir, IMfiles{ImIx,:}(1:end-4),'_',num2str(wluse(i)),'_aligned.tif'],...
                        [IMpaths{ImIx,:}, filesep, IMfiles{ImIx,:}(1:end-4),'_',num2str(wluse(i)),'_aligned.tif']);
                    %         else
                    %             copyfile([pname,filesep,handles.Data.ImageFilename(1:end-4),'_',num2str(wluse(i)),'.tif'],...
                    %                 [pname,filesep,'aligned',filesep,handles.Data.ImageFilename(1:end-4),'_',num2str(wluse(i)),'.tif']);
%                         imwrite(ImAlign{i}(:,:,j),[IMpaths{ImIx,:},filesep,'aligned',filesep,IMfiles{ImIx,:}(1:end-4),'_',num2str(wluse(i)),'.tif'],'WriteMode','append');
%                     end
%                 else
%                     copyfile([IMpaths{ImIx,:},filesep,IMfiles{ImIx,:}(1:end-4),'_',num2str(wluse(i)),'.tif'],...
%                         [IMpaths{ImIx,:},filesep,'aligned',filesep,IMfiles{ImIx,:}(1:end-4),'_',num2str(wluse(i)),'.tif']);
                end
            end
        end
        guidata(hObject,handles);
        
        
        
        cd(Dcur);
        set(handles.StatusText,'String','Batch Processing Complete');
            drawnow
    else
        set(handles.StatusText,'String','Batch Processing Aborted');
        drawnow
        
    end
end


% --- Executes on button press in CutFirstFrame.
function CutFirstFrame_Callback(hObject, eventdata, handles)
% hObject    handle to CutFirstFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CutFirstFrame
