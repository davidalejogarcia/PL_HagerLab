function varargout = JDH1_Analysis_GUI(varargin)
% JDH1_ANALYSIS_GUI M-file for JDH1_Analysis_GUI.fig
%      JDH1_ANALYSIS_GUI, by itself, creates a new JDH1_ANALYSIS_GUI or raises the existing
%      singleton*.
%
%      H = JDH1_ANALYSIS_GUI returns the handle to a new JDH1_ANALYSIS_GUI or the handle to
%      the existing singleton*.
%
%      JDH1_ANALYSIS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in JDH1_ANALYSIS_GUI.M with the given input arguments.
%
%      JDH1_ANALYSIS_GUI('Property','Value',...) creates a new JDH1_ANALYSIS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before JDH1_Analysis_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to JDH1_Analysis_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help JDH1_Analysis_GUI

% Last Modified by GUIDE v2.5 15-Jan-2012 12:17:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @JDH1_Analysis_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @JDH1_Analysis_GUI_OutputFcn, ...
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


% --- Executes just before JDH1_Analysis_GUI is made visible.
function JDH1_Analysis_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to JDH1_Analysis_GUI (see VARARGIN)

% Choose default command line output for JDH1_Analysis_GUI
handles.output = hObject;

Tracks = varargin{1};                    % Import tracks.
Parameters = varargin{2};                % Import stacks of images.
frameTime = varargin{3};                 % Import FrameTime
maxJump = varargin{4};
handles.FileName = varargin{5};
handles.ROIname = varargin{6};

binSize = Parameters.JHistBinSize;
maxFrame = Parameters.JHistMaxFrameN;

binEdges = 0:binSize:maxJump;

% Initial value for the diffusion coefficient
D0 = 1;

% Calculate time-dependent histogram of displacements
[tlist, rlist, JDH] = ...
   calculateJDH(Tracks,maxFrame, binEdges,frameTime, 0, 0);

% Isolate first histogram of displacements;
FJH = [];
FJH(:,1) = rlist;
FJH(:,2) = JDH(1,:);

% Fit first histogram of displacements;
[FJH_Coef1Cmp,FJH_Fit1Cmp,FJH_Sigma1Cmp] = JDHfixedT_1cmp_fit(FJH,frameTime, D0);
[FJH_Coef2Cmp,FJH_Fit2Cmp, FJH_Sigma2Cmp] = JDHfixedT_2cmp_fit(FJH,frameTime, [10*D0, D0]);
[FJH_Coef3Cmp,FJH_Fit3Cmp, FJH_Sigma3Cmp] = JDHfixedT_3cmp_fit(FJH,frameTime, [10*D0, D0, 0.1*D0]);

% Add fits of first histogram of displeacements to FJH.
FJH(:,3) = FJH_Fit1Cmp(:,2);
FJH(:,4) = FJH_Fit2Cmp(:,2);
FJH(:,5) = FJH_Fit3Cmp(:,2);

% Plot the time dependent histogram of displacements in the axes;
axes(handles.axes1);
imagesc(rlist, tlist, JDH);
set(gca, 'FontSize', 10);

xlabel('Jump distance [\mum]');
ylabel('Time [s]');
title('Time-dependent histogram of displacements');
colorbar('location','EastOutside','XTick',0.5,'XTickLabel','Counts');

axis ij

% Pass the output to handles
figname = get(handles.figure1,'Name');
figname = [figname ' ' handles.ROIname];
set(handles.figure1,'Name',figname);

handles.JDH = zeros(maxFrame + 1, length(binEdges));
handles.JDH(2:end,2:end) = JDH;
handles.JDH(2:end, 1) = tlist;
handles.JDH(1, 2:end) = rlist;
handles.FJH = FJH;
handles.FitPar = [FJH_Coef1Cmp, FJH_Coef2Cmp, FJH_Coef3Cmp;...
                   FJH_Sigma1Cmp, 0, FJH_Sigma2Cmp, 0, FJH_Sigma3Cmp, 0];
handles.frameTime =frameTime;
handles.binSize = binSize;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes JDH1_Analysis_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = JDH1_Analysis_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Export Output to main GUI
varargout{1} = handles.JDH;
varargout{2} = handles.FJH;
varargout{3} = handles.FitPar;

% --- Executes on selection change in SelectFit.
function SelectFit_Callback(hObject, eventdata, handles)
% hObject    handle to SelectFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns SelectFit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SelectFit

frameTime = handles.frameTime;
FJH = handles.FJH;
binSize = handles.binSize;
FitPar = handles.FitPar;
JDH = handles.JDH(2:end, 2:end);

tlist = handles.JDH(2:end, 1);
rlist = handles.JDH(1, 2:end)';
rlist(1,2) =tlist(1);




Selection = get(hObject,'Value');

if Selection == 1

    imagesc(rlist(:,1), tlist, JDH);
set(gca, 'FontSize', 10);

xlabel('Jump distance [\mum]');
ylabel('Time [s]');
title('Time-dependent histogram of displacements');
colorbar('location','EastOutside','XTick',0.5,'XTickLabel','Counts');

axis ij
FitString = 'Fit Results';
set(handles.ShowFitResults, 'String', FitString);
    
elseif Selection == 2 % One component diffusion fit

    % plot the histogram of displacements and the fit in the current
    % axes
    
    set(gca, 'FontSize', 10);
    bar(FJH(:,1),FJH(:,2),'EdgeColor', [0 0 0], ...
        'FaceColor','none', 'BarWidth', 1, 'LineWidth', 0.75);
    hold on;
    plot(FJH(:,1),FJH(:,3),'r', 'LineWidth', 1);
    hold off;
    xlabel('Jump distance [\mum]', 'FontSize', 10);
    ylabel('Counts', 'FontSize', 10);
    title({'Fit of the first histogram of displacements';...
        'One Component diffusion Fit';...
        ['\DeltaT = ', num2str(frameTime), ' s']},'FontSize', 10);
    legend({'Displacement Histogram';'One Component Fit'}, 'FontSize', 10);
    xlim ([0 max(FJH(:,1)) + binSize/2]);
    
    % Display fitting results
    
    FitString{1} = 'Estimated Parameters of the one component diffusion fit:';
    FitString{2} = ['n = ' num2str(FitPar(1,1),3),' +/- ', num2str(FitPar(2,1),2)];
    FitString{3} = ['D = ' num2str(FitPar(1,2),3),' +/- ', num2str(FitPar(2,2),2), ' um2/s'];
    FitString{4} = ['SSR = ' num2str(FitPar(1,3));];
    set(handles.ShowFitResults, 'String', FitString);
    
    elseif Selection == 3 % Two component diffusion fit

    % plot the histogram of displacements and the fit in the current
    % axes
    
    FirstCMP = ...
        JDHfixedT_1cmp_fun([FitPar(1,4)* FitPar(1,7), FitPar(1,5)],rlist);
    SecondCMP = ...
        JDHfixedT_1cmp_fun([FitPar(1,4)* (1 - FitPar(1,7)), FitPar(1,6)],rlist);
    
    
    set(gca, 'FontSize', 10);
    bar(FJH(:,1),FJH(:,2),'EdgeColor', [0 0 0], ...
        'FaceColor','none', 'BarWidth', 1, 'LineWidth', 0.75);
    hold on;
    plot(FJH(:,1),FJH(:,4),'r', 'LineWidth', 1);
    plot(FJH(:,1),FirstCMP, '--k','LineWidth', 1);
    plot(FJH(:,1),SecondCMP, '-.k','LineWidth', 1);
    hold off;
    xlabel('Jump distance [\mum]', 'FontSize', 10);
    ylabel('Counts', 'FontSize', 10);
    title({'Fit of the first histogram of displacements';...
        'Two Components diffusion Fit';...
        ['\DeltaT = ', num2str(frameTime), ' s']},'FontSize', 10);
    legend({'Displacement Histogram';'Two Components Fit';...
        'First Component'; 'Second Component'}, 'FontSize', 10);
    xlim ([0 max(FJH(:,1)) + binSize/2]);
    
    % Display fitting results
    
    FitString{1} = 'Estimated Parameters of the two components diffusion fit:';
    FitString{2} = ['n = ' num2str(FitPar(1,4),3),' +/- ', num2str(FitPar(2,4),2)];
    FitString{3} = ['D1 = ' num2str(FitPar(1,5),3),' +/- ', num2str(FitPar(2,5),2), ' um2/s'];
    FitString{4} = ['D2 = ' num2str(FitPar(1,6),3),' +/- ', num2str(FitPar(2,6),2), ' um2/s'];
    FitString{5} = ['f1 = ' num2str(FitPar(1,7),3),' +/- ', num2str(FitPar(2,7),2)];
    FitString{6} = ['SSR = ' num2str(FitPar(1,8));];
    set(handles.ShowFitResults, 'String', FitString);
    
    elseif Selection == 4 % Three component diffusion fit

    % plot the histogram of displacements and the fit in the current
    % axes
    
    FirstCMP = ...
        JDHfixedT_1cmp_fun([FitPar(1,9)* FitPar(1,13), FitPar(1,10)],rlist);
    SecondCMP = ...
        JDHfixedT_1cmp_fun([FitPar(1,9)* FitPar(1,14), FitPar(1,11)],rlist);
    ThirdCMP = JDHfixedT_1cmp_fun...
        ([FitPar(1,9)* (1 - FitPar(1,13) - FitPar(1,14)), FitPar(1,12)],rlist);
    
    
    set(gca, 'FontSize', 10);
    bar(FJH(:,1),FJH(:,2),'EdgeColor', [0 0 0], ...
        'FaceColor','none', 'BarWidth', 1, 'LineWidth', 0.75);
    hold on;
    plot(FJH(:,1),FJH(:,5),'r', 'LineWidth', 1);
    plot(FJH(:,1),FirstCMP, '--k','LineWidth', 1);
    plot(FJH(:,1),SecondCMP, '-.k','LineWidth', 1);
    plot(FJH(:,1),ThirdCMP, ':k','LineWidth', 1);
    hold off;
    xlabel('Jump distance [\mum]', 'FontSize', 10);
    ylabel('Counts', 'FontSize', 10);
    title({'Fit of the first histogram of displacements';...
        'Three Components diffusion Fit';...
        ['\DeltaT = ', num2str(frameTime), ' s']},'FontSize', 10);
    legend({'Displacement Histogram';'Three Components Fit';...
        'First Component'; 'Second Component';'Third Component'}, 'FontSize', 10);
    xlim ([0 max(FJH(:,1)) + binSize/2]);
    
    % Display fitting results
    
    FitString{1} = 'Estimated Parameters of the three components diffusion fit:';
    FitString{2} = ['n = ' num2str(FitPar(1,9),3),' +/- ', num2str(FitPar(2,9),2)];
    FitString{3} = ['D1 = ' num2str(FitPar(1,10),3),' +/- ', num2str(FitPar(2,10),2), ' um2/s'];
    FitString{4} = ['D2 = ' num2str(FitPar(1,11),3),' +/- ', num2str(FitPar(2,11),2), ' um2/s'];
    FitString{5} = ['D3 = ' num2str(FitPar(1,12),3),' +/- ', num2str(FitPar(2,12),2), ' um2/s'];
    FitString{6} = ['f1 = ' num2str(FitPar(1,13),3),' +/- ', num2str(FitPar(2,13),2)];
    FitString{7} = ['f2 = ' num2str(FitPar(1,14),3),' +/- ', num2str(FitPar(2,14),2)];
    FitString{8} = ['SSR = ' num2str(FitPar(1,15));];
    set(handles.ShowFitResults, 'String', FitString);    
end
    


% --- Executes during object creation, after setting all properties.
function SelectFit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelectFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%This button saves the figure currently shown.

frameTime = handles.frameTime;
FJH = handles.FJH;
binSize = handles.binSize;
FitPar = handles.FitPar;
JDH = handles.JDH(2:end, 2:end);

DefaultName = handles.FileName;
SearchStr = '(.*)\.\w*';

tlist = handles.JDH(2:end, 1);
rlist = handles.JDH(1, 2:end)';
rlist(1,2) =tlist(1);

scrsz = get(0,'ScreenSize');
h1 = figure('Position',[.15*scrsz(3) .15*scrsz(4) .4*scrsz(3) .5*scrsz(4)]) ;

Selection = get(handles.SelectFit, 'Value');

if Selection == 1;
    
    imagesc(rlist(:,1), tlist, JDH);
    set(gca, 'FontSize', 10);
    
    xlabel('Jump distance [\mum]');
    ylabel('Time [s]');
    title('Time-dependent histogram of displacements');
    colorbar('location','EastOutside','XTick',0.5,'XTickLabel','Counts');
    
    axis ij
    
    % Save the figures
    
    DefaultName = regexprep(DefaultName, SearchStr, ['$1_' handles.ROIname '_JDH']);
    FilterSpec = {'*.fig';'*.ai';'*.bpm';'*.eps';'*.jpeg';'*.pdf';'*.tif'};
    [FileNameOut,PathNameOut] = uiputfile(FilterSpec,'Save Figure',DefaultName);
    if FileNameOut ~= 0;
        saveas(h1, [PathNameOut,FileNameOut]);
    end
    
elseif Selection == 2 % One component diffusion fit
    
    % plot the histogram of displacements and the fit in the current
    % axes
    
    set(gca, 'FontSize', 10);
    bar(FJH(:,1),FJH(:,2),'EdgeColor', [0 0 0], ...
        'FaceColor','none', 'BarWidth', 1, 'LineWidth', 0.75);
    hold on;
    plot(FJH(:,1),FJH(:,3),'r', 'LineWidth', 1);
    hold off;
    xlabel('Jump distance [\mum]', 'FontSize', 10);
    ylabel('Counts', 'FontSize', 10);
    
    FitString{1} = 'One component diffusion fit of the first Jump Histogram';
    FitString{2} = ['n = ' num2str(FitPar(1,1),3),' \pm ', num2str(FitPar(2,1),2),'; ',...
        'D = ', num2str(FitPar(1,2),3),' \pm ', num2str(FitPar(2,2),2), ' \mum^2/s; '];
    FitString{3}= ['SSR = ' num2str(FitPar(1,3)),'; \DeltaT = ', num2str(frameTime)];
    
    
    title(FitString,'FontSize', 10);
    legend({'Displacement Histogram';'One Component Fit'}, 'FontSize', 10);
    xlim ([0 max(FJH(:,1)) + binSize/2]);
    
    DefaultName = regexprep(DefaultName, SearchStr, ['$1_' handles.ROIname '_FJH_1CMP']);
    FilterSpec = {'*.fig';'*.ai';'*.bpm';'*.eps';'*.jpeg';'*.pdf';'*.tif'};
    [FileNameOut,PathNameOut] = uiputfile(FilterSpec,'Save Figure',DefaultName);
    if FileNameOut ~= 0;
        saveas(h1, [PathNameOut,FileNameOut]);
    end
    
elseif Selection == 3 % Two component diffusion fit
    
    % plot the histogram of displacements and the fit in the current
    % axes
    
    FirstCMP = ...
        JDHfixedT_1cmp_fun([FitPar(1,4)* FitPar(1,7), FitPar(1,5)],rlist);
    SecondCMP = ...
        JDHfixedT_1cmp_fun([FitPar(1,4)* (1 - FitPar(1,7)), FitPar(1,6)],rlist);
    
    
    set(gca, 'FontSize', 10);
    bar(FJH(:,1),FJH(:,2),'EdgeColor', [0 0 0], ...
        'FaceColor','none', 'BarWidth', 1, 'LineWidth', 0.75);
    hold on;
    plot(FJH(:,1),FJH(:,4),'r', 'LineWidth', 1);
    plot(FJH(:,1),FirstCMP, '--k','LineWidth', 1);
    plot(FJH(:,1),SecondCMP, '-.k','LineWidth', 1);
    hold off;
    xlabel('Jump distance [\mum]', 'FontSize', 10);
    ylabel('Counts', 'FontSize', 10);
    FitString{1} = 'Two components diffusion fit of the first Jump Histogram';
    FitString{2} = ['n = ' num2str(FitPar(1,4),3),' \pm ', num2str(FitPar(2,4),2),'; ',...
        'D_1 = ', num2str(FitPar(1,5),3),' \pm ', num2str(FitPar(2,5),2), ' \mum^2/s; ',...
        'D_2 = ', num2str(FitPar(1,6),3),' \pm ', num2str(FitPar(2,6),2), ' \mum^2/s; '];
    FitString{3}= ['f_1 = ', num2str(FitPar(1,7),3),' \pm ', num2str(FitPar(2,7),2), ...
        '; SSR = ' num2str(FitPar(1,8)),'; \DeltaT = ', num2str(frameTime)];
    
    
    title(FitString,'FontSize', 10);
    legend({'Displacement Histogram';'Two Components Fit';...
        'First Component'; 'Second Component'}, 'FontSize', 10);
    xlim ([0 max(FJH(:,1)) + binSize/2]);
    
    DefaultName = regexprep(DefaultName, SearchStr, ['$1_' handles.ROIname '_FJH_2CMP']);
    FilterSpec = {'*.fig';'*.ai';'*.bpm';'*.eps';'*.jpeg';'*.pdf';'*.tif'};
    [FileNameOut,PathNameOut] = uiputfile(FilterSpec,'Save Figure',DefaultName);
    if FileNameOut ~= 0;
        saveas(h1, [PathNameOut,FileNameOut]);
    end
    
elseif Selection == 4 % Three component diffusion fit
    
    % plot the histogram of displacements and the fit in the current
    % axes
    
    FirstCMP = ...
        JDHfixedT_1cmp_fun([FitPar(1,9)* FitPar(1,13), FitPar(1,10)],rlist);
    SecondCMP = ...
        JDHfixedT_1cmp_fun([FitPar(1,9)* FitPar(1,14), FitPar(1,11)],rlist);
    ThirdCMP = JDHfixedT_1cmp_fun...
        ([FitPar(1,9)* (1 - FitPar(1,13) - FitPar(1,14)), FitPar(1,12)],rlist);
    
    
    set(gca, 'FontSize', 10);
    bar(FJH(:,1),FJH(:,2),'EdgeColor', [0 0 0], ...
        'FaceColor','none', 'BarWidth', 1, 'LineWidth', 0.75);
    hold on;
    plot(FJH(:,1),FJH(:,5),'r', 'LineWidth', 1);
    plot(FJH(:,1),FirstCMP, '--k','LineWidth', 1);
    plot(FJH(:,1),SecondCMP, '-.k','LineWidth', 1);
    plot(FJH(:,1),ThirdCMP, ':k','LineWidth', 1);
    hold off;
    xlabel('Jump distance [\mum]', 'FontSize', 10);
    ylabel('Counts', 'FontSize', 10);
    
    FitString{1} = 'Three components diffusion fit of the first Jump Histogram';
    FitString{2} = ['n = ' num2str(FitPar(1,9),3),' \pm ', num2str(FitPar(2,9),2),'; ',...
        'D_1 = ', num2str(FitPar(1,10),3),' \pm ', num2str(FitPar(2,10),2), ' \mum^2/s; ',...
        'D_2 = ', num2str(FitPar(1,11),3),' \pm ', num2str(FitPar(2,11),2), ' \mum^2/s; ',...
        'D_3 = ', num2str(FitPar(1,12),3),' \pm ', num2str(FitPar(2,12),2), ' \mum^2/s; '];
    FitString{3}= ['f_1 = ', num2str(FitPar(1,13),3),' \pm ', num2str(FitPar(2,13),2), ...
        '; f_2 = ', num2str(FitPar(1,14),3),' \pm ', num2str(FitPar(2,14),2),...
        '; SSR = ' num2str(FitPar(1,15)),'; \DeltaT = ', num2str(frameTime)];
    
    title(FitString,'FontSize', 10);
    legend({'Displacement Histogram';'Three Components Fit';...
        'First Component'; 'Second Component';'Third Component'}, 'FontSize', 10);
    xlim ([0 max(FJH(:,1)) + binSize/2]);
    
    DefaultName = regexprep(DefaultName, SearchStr, ['$1_' handles.ROIname '_FJH_3CMP']);
    FilterSpec = {'*.fig';'*.ai';'*.bpm';'*.eps';'*.jpeg';'*.pdf';'*.tif'};
    [FileNameOut,PathNameOut] = uiputfile(FilterSpec,'Save Figure',DefaultName);
    if FileNameOut ~= 0;
        saveas(h1, [PathNameOut,FileNameOut]);
    end
    
    
    
end


% --- Executes on button press in CopyFit.
function CopyFit_Callback(hObject, eventdata, handles)
% hObject    handle to CopyFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


num2clip(handles.FJH);
msgbox({'The Array containing the fitted first jump histogram';...
    'has been copied to the clipboard'});


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

num2clip(handles.FitPar);
msgbox({'The Array containing the Fit Parameters';...
    'has been copied to the clipboard'});

% --- Executes on button press in Copy2DHist.
function Copy2DHist_Callback(hObject, eventdata, handles)
% hObject    handle to Copy2DHist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

num2clip(handles.JDH);
msgbox({'The Array containing the time-dependent histogram';...
    'of displacements has been copied to the clipboard'});


% --- Executes on button press in DoneButton.
function DoneButton_Callback(hObject, eventdata, handles)
% hObject    handle to DoneButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close();

function ShowFitResults_Callback(hObject, eventdata, handles)
% hObject    handle to ShowFitResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ShowFitResults as text
%        str2double(get(hObject,'String')) returns contents of ShowFitResults as a double


% --- Executes during object creation, after setting all properties.
function ShowFitResults_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ShowFitResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
