function varargout = CalcNmin_GUI(varargin)
% CALCNMIN_GUI MATLAB code for CalcNmin_GUI.fig
%      CALCNMIN_GUI, by itself, creates a new CALCNMIN_GUI or raises the existing
%      singleton*.
%
%      H = CALCNMIN_GUI returns the handle to a new CALCNMIN_GUI or the handle to
%      the existing singleton*.
%
%      CALCNMIN_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALCNMIN_GUI.M with the given input arguments.
%
%      CALCNMIN_GUI('Property','Value',...) creates a new CALCNMIN_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CalcNmin_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CalcNmin_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CalcNmin_GUI

% Last Modified by GUIDE v2.5 25-May-2016 13:23:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CalcNmin_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @CalcNmin_GUI_OutputFcn, ...
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


% --- Executes just before CalcNmin_GUI is made visible.
function CalcNmin_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CalcNmin_GUI (see VARARGIN)

% Choose default command line output for CalcNmin_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CalcNmin_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CalcNmin_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function RmaxEdit_Callback(hObject, eventdata, handles)
% hObject    handle to RmaxEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RmaxEdit as text
%        str2double(get(hObject,'String')) returns contents of RmaxEdit as a double


% --- Executes during object creation, after setting all properties.
function RmaxEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RmaxEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DiffCoefEdit_Callback(hObject, eventdata, handles)
% hObject    handle to DiffCoefEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DiffCoefEdit as text
%        str2double(get(hObject,'String')) returns contents of DiffCoefEdit as a double


% --- Executes during object creation, after setting all properties.
function DiffCoefEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DiffCoefEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_intEdit_Callback(hObject, eventdata, handles)
% hObject    handle to T_intEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_intEdit as text
%        str2double(get(hObject,'String')) returns contents of T_intEdit as a double


% --- Executes during object creation, after setting all properties.
function T_intEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_intEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PthEdit_Callback(hObject, eventdata, handles)
% hObject    handle to PthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PthEdit as text
%        str2double(get(hObject,'String')) returns contents of PthEdit as a double


% --- Executes during object creation, after setting all properties.
function PthEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CalcButton.
function CalcButton_Callback(hObject, eventdata, handles)
% hObject    handle to CalcButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rmax = str2double(get(handles.RmaxEdit,'String'));
T_int = str2double(get(handles.T_intEdit,'String'));
D = str2double(get(handles.DiffCoefEdit,'String'));
P_th = str2double(get(handles.PthEdit,'String'));

[Nmin,~] = calculateNmin(rmax,T_int,D,P_th);

set(handles.NminText,'String',num2str(Nmin));
