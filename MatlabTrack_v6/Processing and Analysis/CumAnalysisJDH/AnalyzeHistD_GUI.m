function varargout = AnalyzeHistD_GUI(varargin)
% ANALYZEHISTD_GUI MATLAB code for AnalyzeHistD_GUI.fig
%      ANALYZEHISTD_GUI, by itself, creates a new ANALYZEHISTD_GUI or raises the existing
%      singleton*.
%
%      H = ANALYZEHISTD_GUI returns the handle to a new ANALYZEHISTD_GUI or the handle to
%      the existing singleton*.
%
%      ANALYZEHISTD_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYZEHISTD_GUI.M with the given input arguments.
%
%      ANALYZEHISTD_GUI('Property','Value',...) creates a new ANALYZEHISTD_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AnalyzeHistD_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AnalyzeHistD_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AnalyzeHistD_GUI

% Last Modified by GUIDE v2.5 09-Feb-2012 17:50:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @AnalyzeHistD_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @AnalyzeHistD_GUI_OutputFcn, ...
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


% --- Executes just before AnalyzeHistD_GUI is made visible.
function AnalyzeHistD_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AnalyzeHistD_GUI (see VARARGIN)

% Choose default command line output for AnalyzeHistD_GUI
handles.output = hObject;

%Get input

handles.PathName = varargin{2};
handles.FileNames = varargin{3};
handles.AnalysisS = varargin{4};






% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AnalyzeHistD_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AnalyzeHistD_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



JHistBinSize = str2num(handles.AnalysisS{1});
JHistMaxFrameN = str2num(handles.AnalysisS{2});
maxJump = str2num(handles.AnalysisS{3});
frameTime = str2num(handles.AnalysisS{4});

binEdges = 0:JHistBinSize:maxJump;



% If only one file has been selected change the fileNames variable
% to a cell array

if ~iscell(handles.FileNames)
    handles.FileNames = {handles.FileNames};
end;

nFiles = length(handles.FileNames);

% get the tracks of the selected files
% open each mat-files and retrieve useful information
ROIChoice = {};
for i = 1:nFiles
    
    Temp = load([handles.PathName,handles.FileNames{i}],'Results');
    ROIidx = 0;
    
    if isfield(Temp.Results, 'PreAnalysis') && ...
            isfield(Temp.Results.PreAnalysis, 'Tracks_um');
        Tracks = Temp.Results.PreAnalysis.Tracks_um;
        NParticles{i} = Temp.Results.PreAnalysis.NParticles;
        if isfield(Temp.Results,'Process') %check that the strcutures exist
            if isfield(Temp.Results.Process,'ROIlabel')
                if size(Temp.Results.Process.ROIlabel,1) > 1
                    if ~isempty(ROIChoice) %if a named ROI was selected before check to see if an ROI with that name exists in this dataset
                        if ~strcmp(ROIChoice{1},'All')
                            for j = 1:length(ROIChoice)
                                for k = 1:length(Temp.Results.Process.ROIlabel)
                                    if strcmpi(ROIChoice{j},Temp.Results.Process.ROIlabel{k})
                                        ROIidx = k;
                                        break
                                    end
                                end
                            end
                        else
                            ROIidx = size(Temp.Results.Process.ROIlabel,1) + 1; % If all was selected before, just set the index to something greater than the number of ROIs
                        end
                    end
                    
                    if ROIidx == 0 %if nothing has been found automatically, provide the user with a dialog box to select the ROI
                        
                        ROIstring = Temp.Results.Process.ROIlabel;
                        ROIstring{end+1,1} = 'All';
                        ROIidx = ROIchooseDlg(ROIstring);
                        
                    end
                    %Update the Tracks & NParticles data
                    if ROIidx <= size(Temp.Results.Process.ROIlabel,1)
                        if size(Tracks,2) < 8
                            Tracks = Tracks(Tracks(:,5) == ROIidx,:);
                        else
                            Tracks = Tracks(Tracks(:,9) == ROIidx,:);
                        end
                        NParticles{i} = [NParticles{i}(:,1) NParticles{i}(:,ROIidx+1)];
                        ROIChoice{end+1,1} = Temp.Results.Process.ROIlabel{ROIidx};
                    else
                        tmp = NParticles{i}(:,2:size(Temp.Results.Process.ROIlabel,1)+1);
                        NParticles{i} = [NParticles{i}(:,1) sum(tmp,2)];
                        ROIChoice{end+1,1} = 'All';
                    end
                else
                    ROIChoice{1} = Temp.Results.Process.ROIlabel{1};
                end
            else
                ROIChoice{1} = '';
            end
        end
        [handles.tlist, handles.rlist, JDH(:,:,i)] = ...
            calculateJDH(Tracks,JHistMaxFrameN , binEdges,frameTime, 0, 0);
        
    else
        disp(handles.FileNames{i});
        errordlg('Some of the files have not been preprocessed');
        
        return
    end
end

% Accumulate Jump Histogram
JDH = sum(JDH, 3);

% Calculate the total number of particles;

% find the data containing the lowest number of frames;

nMax = length(NParticles{1}(:,1));
TimePoints = NParticles{1}(:,1);
for i = 1:nFiles
    
    if length(NParticles{i}(:,1)) < nMax;
        
        TimePoints = NParticles{i}(:,1);
        nMax = length(TimePoints);
    end
end

% Accumulate the histogram
Hist_Matrix = zeros(nMax,nFiles);

for i = 1:nFiles
    Hist_Matrix(1:nMax,i) = NParticles{i}(1:nMax,2);
end

CumNParticles(:,1) = TimePoints * frameTime;
CumNParticles(:,2) = sum(Hist_Matrix, 2);

% Fit a double exponential to the photobleaching curve
disp(' ')
disp('____________________________________')
disp('Estimating bleaching characteristics')
disp('____________________________________')

[BleachRates, Dummy, CumNParticles(:,3:6)] =...
    ExpDecay_2Cmp_fit(CumNParticles, [1 0.1]);
disp(['Bleach Rate 1: ', num2str(BleachRates(1), 3), ' s^-1'])
disp(['Bleach Rate 2: ', num2str(BleachRates(2), 3), ' s^-1'])
disp(['Fraction 1: ', num2str(BleachRates(3), 3)])
disp('____________________________________')
disp(' ')
CumNParticles = CumNParticles(:,1:4);
% Plot the photobleaching curve

plot(CumNParticles(:,1), CumNParticles(:,2)/max(CumNParticles(:,4)), '+', 'MarkerEdgeColor', [0.5 0.5 0.5]');
hold on


plot(CumNParticles(:,3), CumNParticles(:,4)/max(CumNParticles(:,4)),'r');

box on;

title({'Photobleaching kinetics:', ['k_{b1} = ', num2str(BleachRates(1), 3),...
    ' s^{-1}, k_{b2} = ', num2str(BleachRates(2), 3), ' s^{-1}, f_1 = ',  ...
    num2str(BleachRates(3), 3)]})

xlabel('Time [s]');
ylabel('Normalized Counts');

hold off


% Ask if you want to correct the jump histogram for photobleaching
answer = questdlg('Do you want to correct the histogram for photobleaching?');

switch answer
    case 'Yes'
        JDH = correctBleachJDH...
            (JDH, BleachRates, handles.tlist, handles.rlist, 0);
end

% Ask if you want to normalize the jump histogram to 1.
answer = questdlg('Do you want to normalize the jump histogram to 1?');

switch answer
    case 'Yes'
        JDH = JDH/sum(JDH(1,:));
end

% Plot heatmap of the histogram of displacements

imagesc(handles.rlist, handles.tlist, JDH);
set(gca, 'FontSize', 10);

xlabel('Jump distance [\mum]');
ylabel('Time [s]');
title('Time-dependent histogram of displacements');
colorbar('location','EastOutside','XTick',0.5,'XTickLabel','Counts');

axis ij

handles.JDH = JDH;

figname = get(gcf,'Name');
figname = [figname ' ' ROIChoice{1}];
set(gcf,'Name',figname);


% Update handles structure
guidata(hObject, handles);
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in MultiD.
function MultiD_Callback(hObject, eventdata, handles)
% hObject    handle to MultiD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in ChooseDraw.
function ChooseDraw_Callback(hObject, eventdata, handles)
% hObject    handle to ChooseDraw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ChooseDraw contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ChooseDraw

JDH = handles.JDH;


Selection =  get(hObject,'Value');

if Selection == 1
    imagesc(handles.rlist, handles.tlist, JDH);
    set(gca, 'FontSize', 10);
    
    xlabel('Jump distance [\mum]');
    ylabel('Time [s]');
    title('Time-dependent histogram of displacements');
    colorbar('location','EastOutside','XTick',0.5,'XTickLabel','Counts');
    
    axis ij
    
elseif Selection == 2
    if isfield(handles, 'fit1Cmp');
    
        Fit = handles.fit1Cmp.Fit;
        FitPar = handles.fit1Cmp.FitPar;
        FitErr = handles.fit1Cmp.FitErr;
        FitSSR = handles.fit1Cmp.SSR;
        
        for t = 1:1:length(handles.tlist);
            tplot = handles.tlist(t)*ones(length(handles.rlist), 1);
            plot3(tplot, handles.rlist, JDH(t,:), 'ok', 'MarkerSize', 4);
            hold on;
            plot3(tplot, handles.rlist, Fit(t,:), 'r');
        end
        hold off;
        
        title({'Jump Histogram Distribution - Full Model Fit',...
            ['D = ', num2str(FitPar(1),3), '\mum^2/s, '],...
            ['k_{on} = ', num2str(FitPar(3),3), 's^{-1} k_{off} = ', num2str(FitPar(4),3), 's^{-1}'],...
            ['dZ =', num2str(FitPar(5),3), '\mum']})
        
        set(gca, 'FontSize', 12);
        view(108,30)
        grid on
        
        xlabel('Time [s]');
        ylabel('Jump Distance [\mum]');
        zlabel('Counts');
    else
        errordlg('You Need to Fit the data first!','Error')
    end
    
elseif Selection == 3
    
    % Plot Fit
    
    if isfield(handles, 'fit2Cmp');
        Fit = handles.fit2Cmp.Fit;
        Fit_2CMP_Ds = handles.fit2Cmp.FitPar(1:2);
        FitPar = handles.fit2Cmp.FitPar(3:end);
        FitErr = handles.fit2Cmp.FitErr;
        FitSSR = handles.fit2Cmp.SSR;
        
        for t = 1:1:length(handles.tlist);
            tplot = handles.tlist(t)*ones(length(handles.rlist), 1);
            plot3(tplot, handles.rlist, JDH(t,:), 'ok', 'MarkerSize', 4);
            hold on;
            plot3(tplot, handles.rlist, Fit(t,:), 'r');
        end
        hold off;
        
        title({'Jump Histogram Distribution - Full Model Fit',...
            ['D_1 = ', num2str(Fit_2CMP_Ds(1)),'\mum^2/s, ', ...
            'D_2 = ', num2str(Fit_2CMP_Ds(2)),'\mum^2/s, '], ...
            ['f_1 = ', num2str(FitPar(1),3), '\mum^2/s, '],...
            ['k_{on} = ', num2str(FitPar(3),3), 's^{-1} k_{off} = ', num2str(FitPar(4),3), 's^{-1}'],...
            ['dZ =', num2str(FitPar(5),3), '\mum']})
        
        
        set(gca, 'FontSize', 12);
        view(108,30)
        grid on
        
        xlabel('Time [s]');
        ylabel('Jump Distance [\mum]');
        zlabel('Counts');
        
    else
        errordlg('You Need to Fit the data first!','Error')
    end
    
    
end







% --- Executes during object creation, after setting all properties.
function ChooseDraw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChooseDraw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OneCMPandbinding.
function OneCMPandbinding_Callback(hObject, eventdata, handles)
% hObject    handle to OneCMPandbinding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the histogram from the handle
JDH = handles.JDH;
rlist = handles.rlist;
tlist = handles.tlist;

% Prompt for the fixed parameters of the fit;
% Get a dialog input to set the parameters;

defaults = {'0.0019'; '0.027'};
prompt = {'Diffusion coefficient of bound molecules [um^2/s]',...
    'Localization accuracy [um]'};

dlgtitle = 'Fixed parameters for the histogram of displacements Fit';
FixedParam = inputdlg(prompt,dlgtitle, 1, defaults);

% Read fixed parameters:
Db = str2num(FixedParam{1});
Sigma = str2num(FixedParam{2});

Par{1} = tlist;
Par{2} = rlist;
Par{3} = Db;
Par{4} = Sigma;

% Prompt for the initial parameters of the fit;

defaults = {'2', '1','.5', '1', '1'};
prompt = {'D_f [um^2/s]',...
    'n', 'k*on [s^-1]','koff [s^-1]', 'dZ [um]'};

dlgtitle = 'Initial guesses for the histogram of displacements Fit';
IniGuess = inputdlg(prompt,dlgtitle, 1, defaults);

Ini(1) = str2num(IniGuess{1});
Ini(2) = str2num(IniGuess{2});
Ini(3) = str2num(IniGuess{3});
Ini(4) = str2num(IniGuess{4});
Ini(5) = str2num(IniGuess{5});

% Fit
[FitPar, FitErr, FitSSR, Fit] = ...
    RD_JD_Fit_dz_BoundDiff(JDH, Par, Ini, 0);



% Plot Fit

for t = 1:1:length(Par{1});
    tplot = Par{1}(t)*ones(length(Par{2}), 1);
    plot3(tplot, Par{2}, JDH(t,:), 'ok', 'MarkerSize', 4);
    hold on;
    plot3(tplot, Par{2}, Fit(t,:), 'r');
end
hold off;

title({'Jump Histogram Distribution - Full Model Fit',...
    ['D = ', num2str(FitPar(1),3), '\mum^2/s, '],...
    ['k_{on} = ', num2str(FitPar(3),3), 's^{-1} k_{off} = ', num2str(FitPar(4),3), 's^{-1}'],...
    ['dZ =', num2str(FitPar(5),3), '\mum']})

set(gca, 'FontSize', 12);
view(108,30)
grid on

xlabel('Time [s]');
ylabel('Jump Distance [\mum]');
zlabel('Counts');




% Prepare output
handles.fit1Cmp.Fit = Fit;
handles.fit1Cmp.FitPar = FitPar;
handles.fit1Cmp.FitErr = FitErr;
handles.fit1Cmp.SSR = FitSSR;

% Set the value of the plot switch.
set(handles.ChooseDraw,'Value', 2);



% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in TwoCMPandbinding.
function TwoCMPandbinding_Callback(hObject, eventdata, handles)
% hObject    handle to TwoCMPandbinding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% Get the histogram from the handle
JDH = handles.JDH;
rlist = handles.rlist;
tlist = handles.tlist;

% Prompt for the fixed parameters of the fit;
% Get a dialog input to set the parameters;

defaults = {'0.0019'; '0.027'};
prompt = {'Diffusion coefficient of bound molecules [um^2/s]',...
    'Localization accuracy [um]'};

dlgtitle = 'Fixed parameters for the histogram of displacements Fit';
FixedParam = inputdlg(prompt,dlgtitle, 1, defaults);

% Read fixed parameters:
Db = str2num(FixedParam{1});
Sigma = str2num(FixedParam{2});

Par{1} = tlist;
Par{2} = rlist;
Par{3} = Db;
Par{4} = Sigma;

% Prompt for the initial parameters of the fit;

defaults = {'5', '.5','.5', '1', '0.5', '1','1'};
prompt = {'D_1 [um^2/s]', 'D_2 [um^2/s]', 'f_1',...
    'n', 'k*on [s^-1]','koff [s^-1]', 'dZ [um]'};

dlgtitle = 'Initial guesses for the histogram of displacements Fit';
IniGuess = inputdlg(prompt,dlgtitle, 1, defaults);


IniD = [str2num(IniGuess{1}), ...
    str2num(IniGuess{2}),str2num(IniGuess{2})/10];
[Fit_2CMP_Ds, Dummy, Fit_2CMP_Ds_Err] = JDHfixedT_3cmp_fit([rlist,JDH(1,:)'] ...
    ,tlist(1), IniD);

Par{5} = Fit_2CMP_Ds(2:3);

Ini(1) = str2num(IniGuess{3});
Ini(2) = str2num(IniGuess{4});
Ini(3) = str2num(IniGuess{5});
Ini(4) = str2num(IniGuess{6});
Ini(5) = str2num(IniGuess{7});

[FitPar, FitErr, FitSSR, Fit] = ...
    RD_JD_Fit_dz_2CMP(JDH, Par, Ini, 0);

% Plot Fit

for t = 1:1:length(Par{1});
    tplot = Par{1}(t)*ones(length(Par{2}), 1);
    plot3(tplot, Par{2}, JDH(t,:), 'ok', 'MarkerSize', 4);
    hold on;
    plot3(tplot, Par{2}, Fit(t,:), 'r');
end
hold off;

title({'Jump Histogram Distribution - Full Model Fit',...
    ['D_1 = ', num2str(Fit_2CMP_Ds(2)),'\mum^2/s, ', ...
    'D_2 = ', num2str(Fit_2CMP_Ds(3)),'\mum^2/s, '], ...
    ['f_1 = ', num2str(FitPar(1),3), '\mum^2/s, '],...
    ['k_{on} = ', num2str(FitPar(3),3), 's^{-1} k_{off} = ', num2str(FitPar(4),3), 's^{-1}'],...
    ['dZ =', num2str(FitPar(5),3), '\mum']})

set(gca, 'FontSize', 12);
view(108,30)
grid on

xlabel('Time [s]');
ylabel('Jump Distance [\mum]');
zlabel('Counts');



% Prepare output
handles.fit2Cmp.Fit = Fit;
handles.fit2Cmp.FitPar = [Fit_2CMP_Ds(2:3), FitPar];
handles.fit2Cmp.FitErr = [Fit_2CMP_Ds_Err(2:3), FitErr];
handles.fit2Cmp.SSR = FitSSR;

% Set the value of the plot switch.
set(handles.ChooseDraw,'Value', 3);


% Update handles structure
guidata(hObject, handles);







% --- Executes on button press in CopyHist.
function CopyHist_Callback(hObject, eventdata, handles)


OUT = zeros(size(handles.JDH) + 1);

OUT(2:end,1) = handles.tlist;
OUT(1, 2:end) = handles.rlist;
OUT(2:end, 2:end) = handles.JDH;

num2clip(OUT);
msgbox('The array with the histogram of displacements has been copied to the clipboard');




% hObject    handle to CopyHist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in CopyFit.
function CopyFit_Callback(hObject, eventdata, handles)

Selection = get(handles.ChooseDraw,'Value');
if Selection == 2
OUT= zeros(size(handles.JDH) + 1);
OUT(2:end,1) = handles.tlist;
OUT(1, 2:end) = handles.rlist;
OUT(2:end, 2:end) = handles.fit1Cmp.Fit;
num2clip(OUT);

msgbox('The array with the 1CMP Full model has been copied to the clipboard');

elseif Selection ==3
OUT= zeros(size(handles.JDH) + 1);
OUT(2:end,1) = handles.tlist;
OUT(1, 2:end) = handles.rlist;
OUT(2:end, 2:end) = handles.fit2Cmp.Fit;
num2clip(OUT);
msgbox('The array with the 2CMP Full model has been copied to the clipboard');


else
    errordlg('Please select which fit do you want to copy');
end





% hObject    handle to CopyFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in CopyFitPar.
function CopyFitPar_Callback(hObject, eventdata, handles)

Selection = get(handles.ChooseDraw,'Value');
if Selection == 2
    OUT = [handles.fit1Cmp.FitPar, handles.fit1Cmp.SSR; ...
        handles.fit1Cmp.FitErr, 0];
        
num2clip(OUT);

msgbox('The array with the 1CMP Full model fit parameters has been copied to the clipboard');

elseif Selection ==3
    

    OUT = [handles.fit2Cmp.FitPar, handles.fit2Cmp.SSR; ...
        handles.fit2Cmp.FitErr, 0];
        
num2clip(OUT);

msgbox('The array with the 2CMP Full model fit parameters has been copied to the clipboard');



else
    errordlg('Please select which fit do you want to copy');
end



% hObject    handle to CopyFitPar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in SaveFigure.
function SaveFigure_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
