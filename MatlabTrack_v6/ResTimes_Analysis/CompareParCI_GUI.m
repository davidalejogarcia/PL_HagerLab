function varargout = CompareParCI_GUI(varargin)
% COMPAREPARCI_GUI MATLAB code for CompareParCI_GUI.fig
%      COMPAREPARCI_GUI, by itself, creates a new COMPAREPARCI_GUI or raises the existing
%      singleton*.
%
%      H = COMPAREPARCI_GUI returns the handle to a new COMPAREPARCI_GUI or the handle to
%      the existing singleton*.
%
%      COMPAREPARCI_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COMPAREPARCI_GUI.M with the given input arguments.
%
%      COMPAREPARCI_GUI('Property','Value',...) creates a new COMPAREPARCI_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CompareParCI_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CompareParCI_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CompareParCI_GUI

% Last Modified by GUIDE v2.5 22-Sep-2017 11:40:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @CompareParCI_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @CompareParCI_GUI_OutputFcn, ...
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


% --- Executes just before CompareParCI_GUI is made visible.
function CompareParCI_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CompareParCI_GUI (see VARARGIN)

% Choose default command line output for CompareParCI_GUI
handles.output = hObject;

handles.Data = [];
handles.Analysis = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CompareParCI_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CompareParCI_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_fname1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fname1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fname1 as text
%        str2double(get(hObject,'String')) returns contents of edit_fname1 as a double

fname = get(hObject,'String');
if ~isempty(fname)
    try f1 = load(fname);
    catch
        errordlg('File not supported')
        return
    end
    [handles.Data.path1,handles.Data.file1] = fileparts(fname);
    handles.Data.file1 = [handles.Data.file1,'.mat'];
    
    handles.Data.fitpar1 = f1.Results.FitPar.Surv_PB;
    
    var1_select = get(handles.popup_var1,'Value');
    if var1_select == 5 || var1_select == 6
        fast_frac = handles.Data.fitpar1(1,5);
        bound_frac = handles.Data.fitpar1(1,6);
        fast_err = handles.Data.fitpar1(2,5);
        bound_err = handles.Data.fitpar1(2,6);
        
        if var1_select == 5
            var1(1,1) = fast_frac*bound_frac;
        else
            var1(1,1) = bound_frac(1 - fast_frac);
        end
        var1(2,1) = (fast_frac*bound_frac)*sqrt((fast_err/fast_frac).^2 +...
            (bound_err/bound_frac).^2);
        str3 = '';
        
        
        
        
    elseif var1_select == 7
        var1 = handles.Data.fitpar1(:,6);
        str3 = '';
    elseif var1_select == 2
        var1 = handles.Data.fitpar1(:,2);
        str3 = '';
    elseif var1_select == 1 || var1_select == 3 || var1_select == 4
        var1_tmp = handles.Data.fitpar1(:,var1_select);
        
        
        var1(1,1) = 1/var1_tmp(1);
        var1(2,1) = (var1(1)/var1_tmp(1))*var1_tmp(2);
        str3 = 's';
    end
    
    
    
    handles.Analysis.Var1_val = var1;
    
    var1_2 = [0;0];
    
    [h0,pval0] = compareCI(var1,var1_2);
    
    str1 = ['Value: ', num2str(var1(1),3), ' +/- ', num2str(var1(2),3), ' ',str3];
    str2 = ['Significantly >0: ',num2str(h0), ', p-value: ', num2str(pval0,3)];
    
    
    set(handles.text_Par1,'String',{str1;str2});
    set(handles.popup_var1,'Enable','on');
    set(handles.popup_var2,'Enable','on');
    set(handles.edit_fname2,'Enable','on');
    set(handles.button_fileOpen2,'Enable','on');
    
    if ~isfield(handles.Data,'fitpar2')
        var2_select = get(handles.popup_var2,'Value');
        if var2_select == 5 || var2_select == 6
            fast_frac = handles.Data.fitpar1(1,5)*handles.Data.fitpar1(1,6);
            diff_frac = handles.Data.fitpar1(1,6);
            fast_err = handles.Data.fitpar1(2,5);
            diff_err = handles.Data.fitpar1(2,6);
            
            if var2_select == 5
                var2(1,1) = fast_frac;
            else
                var2(1,1) = diff_frac - fast_frac;
            end
            var2(2,1) = fast_frac*sqrt((fast_err/fast_frac).^2 +...
                (diff_err/diff_frac).^2);
            str3 = '';
            
            
            
            
        elseif var2_select == 7
            var2 = handles.Data.fitpar1(:,6);
            str3 = '';
        elseif var2_select == 2
            var2 = handles.Data.fitpar1(:,2);
            str3 = '';
        elseif var2_select == 1 || var2_select == 3 || var2_select == 4
            var2_tmp = handles.Data.fitpar1(:,var2_select);
            
            
            var2(1,1) = 1/var2_tmp(1);
            var2(2,1) = (var2(1)/var2_tmp(1))*var2_tmp(2);
            str3 = 's';
        end
        
        
        
        handles.Analysis.Var2_val = var1;
        
        var2_2 = [0;0];
        
        [h0,pval0] = compareCI(var2,var2_2);
        
        str1 = ['Value: ', num2str(var2(1),3), ' +/- ', num2str(var2(2),3), ' ',str3];
        str2 = ['Significantly >0: ',num2str(h0), ', p-value: ', num2str(pval0,3)];
        
        set(handles.text_Par2,'String',{str1;str2});
    end
    set(handles.button_compare,'Enable','on');
end

guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function edit_fname1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fname1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_fileOpen1.
function button_fileOpen1_Callback(hObject, eventdata, handles)
% hObject    handle to button_fileOpen1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles.Data,'path1')
    defPath = handles.Data.path1;
elseif isfield(handles.Data,'path2')
    defPath = handles.Data.path2;
else
    defPath = pwd;
end


[fname,pname] = uigetfile('*.mat','Open Residence Time Analysis Files',defPath);

if fname ~= 0
    try f1 = load(fullfile(pname,fname));
    catch
        errordlg('File not supported')
        return
    end
    handles.Data.path1 = pname;
    handles.Data.file1 = fname;
    
    set(handles.edit_fname1,'String',fullfile(pname,fname));
    
    handles.Data.fitpar1 = f1.Results.FitPar.Surv_PB;
    
    var1_select = get(handles.popup_var1,'Value');
    if var1_select == 5 || var1_select == 6
        fast_frac = handles.Data.fitpar1(1,5);
        bound_frac = handles.Data.fitpar1(1,6);
        fast_err = handles.Data.fitpar1(2,5);
        bound_err = handles.Data.fitpar1(2,6);
        
        if var1_select == 5
            var1(1,1) = fast_frac*bound_frac;
        else
            var1(1,1) = bound_frac(1 - fast_frac);
        end
        var1(2,1) = (fast_frac*bound_frac)*sqrt((fast_err/fast_frac).^2 +...
            (bound_err/bound_frac).^2);
        str3 = '';
        
        
        
        
    elseif var1_select == 7
        var1 = handles.Data.fitpar1(:,6);
        str3 = '';
    elseif var1_select == 2
        var1 = handles.Data.fitpar1(:,2);
        str3 = '';
    elseif var1_select == 1 || var1_select == 3 || var1_select == 4
        var1_tmp = handles.Data.fitpar1(:,var1_select);
        
        
        var1(1,1) = 1/var1_tmp(1);
        var1(2,1) = (var1(1)/var1_tmp(1))*var1_tmp(2);
        str3 = 's';
    end
    
    
    
    handles.Analysis.Var1_val = var1;
    
    var1_2 = [0;0];
    
    [h0,pval0] = compareCI(var1,var1_2);
    
    str1 = ['Value: ', num2str(var1(1),3), ' +/- ', num2str(var1(2),3), ' ',str3];
    str2 = ['Significantly >0: ',num2str(h0), ', p-value: ', num2str(pval0,3)];
    
    
    set(handles.text_Par1,'String',{str1;str2});
    set(handles.popup_var1,'Enable','on');
    set(handles.popup_var2,'Enable','on');
    set(handles.edit_fname2,'Enable','on');
    set(handles.button_fileOpen2,'Enable','on');
    set(handles.menu_copyVar,'Enable','on');
    
    if ~isfield(handles.Data,'fitpar2')
        var2_select = get(handles.popup_var2,'Value');
        if var2_select == 5 || var2_select == 6
            fast_frac = handles.Data.fitpar1(1,5)*handles.Data.fitpar1(1,6);
            diff_frac = handles.Data.fitpar1(1,6);
            fast_err = handles.Data.fitpar1(2,5);
            diff_err = handles.Data.fitpar1(2,6);
            
            if var2_select == 5
                var2(1,1) = fast_frac;
            else
                var2(1,1) = diff_frac - fast_frac;
            end
            var2(2,1) = fast_frac*sqrt((fast_err/fast_frac).^2 +...
                (diff_err/diff_frac).^2);
            str3 = '';
            
            
            
            
        elseif var2_select == 7
            var2 = handles.Data.fitpar1(:,6);
            str3 = '';
        elseif var2_select == 2
            var2 = handles.Data.fitpar1(:,2);
            str3 = '';
        elseif var2_select == 1 || var2_select == 3 || var2_select == 4
            var2_tmp = handles.Data.fitpar1(:,var2_select);
            
            
            var2(1,1) = 1/var2_tmp(1);
            var2(2,1) = (var2(1)/var2_tmp(1))*var2_tmp(2);
            str3 = 's';
        end
        
        
        
        handles.Analysis.Var2_val = var1;
        
        var2_2 = [0;0];
        
        [h0,pval0] = compareCI(var2,var2_2);
        
        str1 = ['Value: ', num2str(var2(1),3), ' +/- ', num2str(var2(2),3), ' ',str3];
        str2 = ['Significantly >0: ',num2str(h0), ', p-value: ', num2str(pval0,3)];
        
        set(handles.text_Par2,'String',{str1;str2});
    end
    set(handles.button_compare,'Enable','on');
end

guidata(hObject,handles);

function edit_fname2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fname2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fname2 as text
%        str2double(get(hObject,'String')) returns contents of edit_fname2 as a double
fname = get(hObject,'String');
if ~isempty(fname)
    try f2 = load(fullfile(pname,fname));
    catch
        errordlg('File not supported')
        return
    end
    [handles.Data.path2,handles.Data.file2] = fileparts(fname);
    handles.Data.file2 = [handles.Data.file2,'.mat'];
    
    handles.Data.fitpar2 = f2.Results.FitPar.Surv_PB;
    
    var2_select = get(handles.popup_var2,'Value');
    if var2_select == 5 || var2_select == 6
        fast_frac = handles.Data.fitpar1(1,5);
        bound_frac = handles.Data.fitpar1(1,6);
        fast_err = handles.Data.fitpar1(2,5);
        bound_err = handles.Data.fitpar1(2,6);
        
        if var2_select == 5
            var2(1,1) = fast_frac*bound_frac;
        else
            var2(1,1) = bound_frac(1 - fast_frac);
        end
        var2(2,1) = (fast_frac*bound_frac)*sqrt((fast_err/fast_frac).^2 +...
            (bound_err/bound_frac).^2);
        str3 = '';
        
        
        
        
    elseif var2_select == 7
        var2 = handles.Data.fitpar2(:,6);
        str3 = '';
    elseif var2_select == 2
        var2 = handles.Data.fitpar2(:,2);
        str3 = '';
    elseif var2_select == 1 || var2_select == 3 || var2_select == 4
        var1_tmp = handles.Data.fitpar2(:,var2_select);
        
        
        var2(1,1) = 1/var1_tmp(1);
        var2(2,1) = (var2(1)/var1_tmp(1))*var1_tmp(2);
        str3 = 's';
    end
    
    
    
    handles.Analysis.Var2_val = var2;
    
    var2_2 = [0;0];
    
    [h0,pval0] = compareCI(var2,var2_2);
    
    str1 = ['Value: ', num2str(var2(1),3), ' +/- ', num2str(var2(2),3), ' ',str3];
    str2 = ['Significantly >0: ',num2str(h0), ', p-value: ', num2str(pval0,3)];
    
    set(handles.text_Par2,'String',{str1;str2});
    
    set(handles.text_Par2,'String',{str1;str2});
    
    
end

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_fname2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fname2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_fileOpen2.
function button_fileOpen2_Callback(hObject, eventdata, handles)
% hObject    handle to button_fileOpen2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles.Data,'path1')
    defPath = handles.Data.path1;
elseif isfield(handles.Data,'path2')
    defPath = handles.Data.path2;
else
    defPath = pwd;
end


[fname,pname] = uigetfile('*.mat','Open Residence Time Analysis Files',defPath);

if fname ~= 0
    try f2 = load(fullfile(pname,fname));
    catch
        errordlg('File not supported')
        return
    end
    handles.Data.path2 = pname;
    handles.Data.file2 = fname;
    
    set(handles.edit_fname2,'String',fullfile(pname,fname));
    
    handles.Data.fitpar2 = f2.Results.FitPar.Surv_PB;
    
    var2_select = get(handles.popup_var2,'Value');
    if var2_select == 5 || var2_select == 6
        fast_frac = handles.Data.fitpar1(1,5);
        bound_frac = handles.Data.fitpar1(1,6);
        fast_err = handles.Data.fitpar1(2,5);
        bound_err = handles.Data.fitpar1(2,6);
        
        if var2_select == 5
            var2(1,1) = fast_frac*bound_frac;
        else
            var2(1,1) = bound_frac(1 - fast_frac);
        end
        var2(2,1) = (fast_frac*bound_frac)*sqrt((fast_err/fast_frac).^2 +...
            (bound_err/bound_frac).^2);
        str3 = '';
        
        
        
        
    elseif var2_select == 7
        var2 = handles.Data.fitpar2(:,6);
        str3 = '';
    elseif var2_select == 2
        var2 = handles.Data.fitpar2(:,2);
        str3 = '';
    elseif var2_select == 1 || var2_select == 3 || var2_select == 4
        var1_tmp = handles.Data.fitpar2(:,var2_select);
        
        
        var2(1,1) = 1/var1_tmp(1);
        var2(2,1) = (var2(1)/var1_tmp(1))*var1_tmp(2);
        str3 = 's';
    end
    
    
    
    handles.Analysis.Var2_val = var2;
    
    var2_2 = [0;0];
    
    [h0,pval0] = compareCI(var2,var2_2);
    
    str1 = ['Value: ', num2str(var2(1),3), ' +/- ', num2str(var2(2),3), ' ',str3];
    str2 = ['Significantly >0: ',num2str(h0), ', p-value: ', num2str(pval0,3)];
    
    set(handles.text_Par2,'String',{str1;str2});
    
    set(handles.text_Par2,'String',{str1;str2});
    
    
end

guidata(hObject,handles);

% --- Executes on selection change in popup_var1.
function popup_var1_Callback(hObject, eventdata, handles)
% hObject    handle to popup_var1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_var1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_var1

var1_select = get(hObject,'Value');
if var1_select == 5 || var1_select == 6
    fast_frac = handles.Data.fitpar1(1,5);
    bound_frac = handles.Data.fitpar1(1,6);
    fast_err = handles.Data.fitpar1(2,5);
    bound_err = handles.Data.fitpar1(2,6);
    
    if var1_select == 5
        var1(1,1) = fast_frac*bound_frac;
    else
        var1(1,1) = bound_frac*(1 - fast_frac);
    end
    var1(2,1) = (fast_frac*bound_frac)*sqrt((fast_err/fast_frac).^2 +...
        (bound_err/bound_frac).^2);
    str3 = '';
    
    
    
    
elseif var1_select == 7
    var1 = handles.Data.fitpar1(:,6);
    str3 = '';
elseif var1_select == 2
    var1 = handles.Data.fitpar1(:,2);
    str3 = '';
elseif var1_select == 1 || var1_select == 3 || var1_select == 4
    var1_tmp = handles.Data.fitpar1(:,var1_select);
    
    
    var1(1,1) = 1/var1_tmp(1);
    var1(2,1) = (var1(1)/var1_tmp(1))*var1_tmp(2);
    str3 = 's';
end



handles.Analysis.Var1_val = var1;

var2 = [0;0];

[h0,pval0] = compareCI(var1,var2);

str1 = ['Value: ', num2str(var1(1),3), ' +/- ', num2str(var1(2),3), ' ',str3];
str2 = ['Significantly >0: ',num2str(h0), ', p-value: ', num2str(pval0,3)];

set(handles.text_Par1,'String',{str1;str2});


guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function popup_var1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_var1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_var2.
function popup_var2_Callback(hObject, eventdata, handles)
% hObject    handle to popup_var2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_var2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_var2

var1_select = get(hObject,'Value');
if isfield(handles.Data,'fitpar2')
    fitpar2 = handles.Data.fitpar2;
else
    fitpar2 = handles.Data.fitpar1;
end

if var1_select == 5 || var1_select == 6
    fast_frac = fitpar2(1,5);
    bound_frac = fitpar2(1,6);
    fast_err = fitpar2(2,5);
    bound_err = fitpar2(2,6);
    
    if var1_select == 5
        var1(1,1) = fast_frac*bound_frac;
    else
        var1(1,1) = bound_frac*(1 - fast_frac);
    end
    var1(2,1) = (fast_frac*bound_frac)*sqrt((fast_err/fast_frac).^2 +...
        (bound_err/bound_frac).^2);
    str3 = '';
    
    
    
    
elseif var1_select == 7
    var1 = fitpar2(:,6);
    str3 = '';
elseif var1_select == 2
    var1 = fitpar2(:,2);
    str3 = '';
elseif var1_select == 1 || var1_select == 3 || var1_select == 4
    var1_tmp = fitpar2(:,var1_select);
    
    
    var1(1,1) = 1/var1_tmp(1);
    var1(2,1) = (var1(1)/var1_tmp(1))*var1_tmp(2);
    str3 = 's';
end



handles.Analysis.Var2_val = var1;

var2 = [0;0];

[h0,pval0] = compareCI(var1,var2);

str1 = ['Value: ', num2str(var1(1),3), ' +/- ', num2str(var1(2),3), ' ',str3];
str2 = ['Significantly >0: ',num2str(h0), ', p-value: ', num2str(pval0,3)];

set(handles.text_Par2,'String',{str1;str2});


guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popup_var2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_var2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_compare.
function button_compare_Callback(hObject, eventdata, handles)
% hObject    handle to button_compare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val1 = handles.Analysis.Var1_val;
val2 = handles.Analysis.Var2_val;


[h,p,D,CI_tot] = compareCI(val1,val2);

str1 = ['Difference: ', num2str(D,3), ' +/- ', num2str(CI_tot,3), ' '];
str2 = ['Significant Difference: ',num2str(h), ', p-value: ', num2str(p,3)];

handles.Analysis.Var1Var2_Diff = [D; CI_tot];
handles.Analysis.Var1Var2_hypTest = [h; p];


set(handles.text_Par1Par2,'String',{str1;str2});
guidata(hObject,handles);


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_copyVar_Callback(hObject, eventdata, handles)
% hObject    handle to menu_copyVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varOut = handles.Analysis.Var1_val;

if isfield(handles.Analysis,'Var2_val')
    varOut = [varOut,handles.Analysis.Var2_val];
end
if isfield(handles.Analysis,'Var1Var2_Diff')
    varOut = [varOut,handles.Analysis.Var1Var2_Diff,handles.Analysis.Var1Var2_hypTest];
end

num2clip(varOut);
msgbox('Values have been copied to the clipboard');
