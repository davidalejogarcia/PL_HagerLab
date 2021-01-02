function varargout = SetROIParamGUI(varargin)
% SETROIPARAMGUI MATLAB code for SetROIParamGUI.fig
%      SETROIPARAMGUI, by itself, creates a new SETROIPARAMGUI or raises the existing
%      singleton*.
%
%      H = SETROIPARAMGUI returns the handle to a new SETROIPARAMGUI or the handle to
%      the existing singleton*.
%
%      SETROIPARAMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SETROIPARAMGUI.M with the given input arguments.
%
%      SETROIPARAMGUI('Property','Value',...) creates a new SETROIPARAMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SetROIParamGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SetROIParamGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SetROIParamGUI

% Last Modified by GUIDE v2.5 14-Dec-2016 09:48:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SetROIParamGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SetROIParamGUI_OutputFcn, ...
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


% --- Executes just before SetROIParamGUI is made visible.
function SetROIParamGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SetROIParamGUI (see VARARGIN)

% Choose default command line output for SetROIParamGUI
handles.output = hObject;

if ~isempty(varargin)
    ROItype = varargin{1};
    ROIsize = varargin{2};
    
    set(handles.ShapePop,'Value',ROItype);
    set(handles.dimValEdit,'String',num2str(ROIsize));
    
    if get(handles.ShapePop,'Value') == 1
        set(handles.dimText,'String','Diameter');
        Area = pi*((ROIsize/2)^2);
    else
        set(handles.dimText,'String','Side');
        Area = (ROIsize^2);
    end
    set(handles.AreaText,'String',num2str(Area,5));
else
    ROItype = 1;
    ROIsize = 14;
end


setappdata(0,'ROItype',ROItype);
setappdata(0,'ROIsize',ROIsize);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SetROIParamGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SetROIParamGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in ShapePop.
function ShapePop_Callback(hObject, eventdata, handles)
% hObject    handle to ShapePop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ShapePop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ShapePop
shapeVal = get(hObject,'Value');
d = str2double(get(handles.dimValEdit,'String'));

if shapeVal == 1
    set(handles.dimText,'String','Diameter');
    Area = pi*((d/2).^2);
else
    set(handles.dimText,'String','Side');
    Area = d.^2;
end
set(handles.AreaText,'String',num2str(Area,5));


setappdata(0,'ROItype',shapeVal);
guidata(hObject,handles);    

% --- Executes during object creation, after setting all properties.
function ShapePop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ShapePop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dimValEdit_Callback(hObject, eventdata, handles)
% hObject    handle to dimValEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dimValEdit as text
%        str2double(get(hObject,'String')) returns contents of dimValEdit as a double
d = str2double(get(hObject,'String'));
shapeVal = get(handles.ShapePop,'Value');

if shapeVal == 1
    Area = pi*((d/2).^2);
else
    
    Area = d.^2;
end
set(handles.AreaText,'String',num2str(Area,5));
setappdata(0,'ROIsize',str2double(get(hObject,'String')));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function dimValEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dimValEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in DoneButton.
function DoneButton_Callback(hObject, eventdata, handles)
% hObject    handle to DoneButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close();



function AreaText_Callback(hObject, eventdata, handles)
% hObject    handle to AreaText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AreaText as text
%        str2double(get(hObject,'String')) returns contents of AreaText as a double
shapeVal = get(handles.ShapePop,'Value');
Area = str2double(get(hObject,'String'));

if shapeVal == 1
    d = 2*sqrt(Area/pi);
else
    
    d = sqrt(Area);
end
set(handles.dimValEdit,'String',num2str(d,3));
setappdata(0,'ROIsize',d);
guidata(hObject,handles);
