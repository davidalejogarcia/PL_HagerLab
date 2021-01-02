function varargout = DiffCoeffMerge(varargin)
% DIFFCOEFFMERGE MATLAB code for DiffCoeffMerge.fig
%      DIFFCOEFFMERGE, by itself, creates a new DIFFCOEFFMERGE or raises the existing
%      singleton*.
%
%      H = DIFFCOEFFMERGE returns the handle to a new DIFFCOEFFMERGE or the handle to
%      the existing singleton*.
%
%      DIFFCOEFFMERGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIFFCOEFFMERGE.M with the given input arguments.
%
%      DIFFCOEFFMERGE('Property','Value',...) creates a new DIFFCOEFFMERGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DiffCoeffMerge_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DiffCoeffMerge_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DiffCoeffMerge

% Last Modified by GUIDE v2.5 15-Sep-2016 17:10:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DiffCoeffMerge_OpeningFcn, ...
                   'gui_OutputFcn',  @DiffCoeffMerge_OutputFcn, ...
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


% --- Executes just before DiffCoeffMerge is made visible.
function DiffCoeffMerge_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DiffCoeffMerge (see VARARGIN)

fileNames = varargin{1};
pathName = varargin{2};

frameTime = varargin{3};
nbins = varargin{4};
handles.frameTime = frameTime;
handles.FileName = fileNames{1};
handles.PathName = pathName;
%initialize vector to store D
D = [];

ROIChoice = {};
ROIorClass = -1; % -1: not chosen, 0: ROI name, 1: Class Name
for i = 1:length(fileNames)
    ROIidx = 0;
    Temp = load(fullfile(pathName,fileNames{i}));
    if isfield(Temp.Results,'Process') %check that the strcutures exist
        if isfield(Temp.Results.Process,'ROIlabel')
            if isfield(Temp.Results.Process,'AllROIClasses') && ROIorClass == -1
                ROIClass =  questdlg('Do you want to separate data based on ROI names or Class Names?','ROI or Class','ROI','Class','Class');
                if strcmp(ROIClass,'ROI');
                    ROIorClass = 0;
                else
                    ROIorClass = 1;
                end
            elseif ~isfield(Temp.Results.Process,'AllROIClasses')
                ROIorClass = 0;
            end
            if ROIorClass == 0
                roiLabels = Temp.Results.Process.ROIlabel;
            elseif ROIorClass == 1
                if ~isfield(Temp.Results.Process,'AllROIClasses')
                    errordlg(['File does not contain Class data: ,' handles.FileNames{i}]);
                else
                    roiLabels = Temp.Results.Process.AllROIClasses;
                end
            end
            if size(roiLabels,1) > 1
                if ~isempty(ROIChoice) %if a named ROI was selected before check to see if an ROI with that name exists in this dataset
                    if ~strcmp(ROIChoice{1},'All')
                        for j = 1:length(ROIChoice)
                            for k = 1:length(roiLabels)
                                if strcmpi(ROIChoice{j},roiLabels{k})
                                    ROIidx = k;
                                    break
                                end
                            end
                        end
                    else
                        ROIidx = size(roiLabels,1) + 1; % If all was selected before, just set the index to something greater than the number of ROIs
                    end
                end
                
                if ROIidx == 0 %if nothing has been found automatically, provide the user with a dialog box to select the ROI
                    
                    ROIstring = roiLabels;
                    ROIstring{end+1,1} = 'All';
                    ROIidx = ROIchooseDlg(ROIstring);
                    
                end
                %Update the Tracks & NParticles data
                if ROIidx <= size(roiLabels,1)
                    if ROIorClass == 0
                        ROIvec = ROIidx;
                    elseif ROIorClass == 1
                        ROIvec = [];
                        for nROI = 1:size(Temp.Results.Process.ROIpos)
                            if strcmpi(Temp.Results.Process.ROIClass{nROI},ROIChoice)
                                ROIvec = [ROIvec;nROI];
                            end
                        end
                    end
                    D_tmp = diffusionImageCalc(Temp.Results,0.8,5,frameTime,ROIvec);
                else
                    D_tmp = diffusionImageCalc(Temp.Results,0.8,5,frameTime);
                    ROIChoice{end+1,1} = 'All';
                end
            else
               ROIChoice{1} = roiLabels{1};
               D_tmp = diffusionImageCalc(Temp.Results,0.8,5,frameTime);
            end
        else
            ROIChoice{1} = '';
        end
    else
        ROIChoice{1} = '';
    end
                    
    
    D = [D; D_tmp];
end
Dlog = log10(D);

%Generate the histogram with at least 25 bins
if nbins == 0
    nbins = max(25,length(Dlog)/8);
end

[n,x] = hist(Dlog,nbins);
handles.binSize = x(2) - x(1);

DLog_hist = [x',n'];

axes(handles.axes1);
bar(DLog_hist(:,1),DLog_hist(:,2),'EdgeColor', [0 0 0], ...
'FaceColor','none', 'BarWidth', 1, 'LineWidth', 0.75);
xlabel('LogD (\mum^2/s)')
ylabel('Frequency');
hold on;
%calculate SST
SST = sum((DLog_hist(:,2) - mean(DLog_hist(:,2))).^2);
%Perform fitting
fit2Show = 0;
handles.Data = DLog_hist;

%1-Comp
[G_Coef1, G_Sigma1, G_Fit1] = Gauss_1Cmp_fit(DLog_hist, [-0.1,0.1]);
FitPar1 = [G_Coef1;G_Sigma1];
FitPar1(1,:) = 10.^FitPar1(1,:);
FitPar1(2,:) = FitPar1(1,:) - 10.^(G_Coef1 - G_Sigma1);

FitHist1 = G_Fit1;

SSR1 = sum((DLog_hist(:,2) - G_Fit1(:,2)).^2);

R2_1 = 1 - (SSR1/SST);
FitPar1(:,3) = [SSR1;R2_1];
if R2_1 >= 0.95
    fit2Show = 1;
end

%2-Comp
k0_2 = [G_Coef1(1:2),G_Coef1(1)/2,G_Coef1(2)];
[G_Coef2, G_Sigma2, G_Fit2] = Gauss_2Cmp_fit(DLog_hist, k0_2);

% %keep the fitted D.F.s in increasing order.
% if G_Coef2(1) > G_Coef2(3)
%     tmp = G_Coef2(1:2);
%     G_Coef2(1:2) = G_Coef2(3:4);
%     G_Coef2(3:4) = tmp;
%     G_Coef2(5) = 1 - G_Coef2(5);
%     tmp = G_Sigma2(1:2);
%     G_Sigma2(1:2) = G_Sigma2(3:4);
%     G_Sigma2(3:4) = tmp;
%     tmp = G_Fit2(:,3);
%     G_Fit2(:,3) = G_Fit2(:,4);
%     G_Fit2(:,4) = tmp;
% end

FitPar2 = [G_Coef2;G_Sigma2];
FitPar2(1,:) = [10.^FitPar2(1,1:4),FitPar2(1,5)];
FitPar2(2,:) = [FitPar2(1,1:4) - 10.^(G_Coef2(1:4) - G_Sigma2(1:4)),FitPar2(2,5)];

FitHist2 = G_Fit2;

SSR2 = sum((DLog_hist(:,2) - G_Fit2(:,2)).^2);

R2_2 = 1 - (SSR2/SST);
FitPar2(:,6) = [SSR2;R2_2];
if R2_2 >= 0.95 && fit2Show == 0
    fit2Show = 2;
end

%3-Comp
k0_3 = [G_Coef2(1:4),G_Coef2(1)/2,G_Coef2(2)];
[G_Coef3, G_Sigma3, G_Fit3] = Gauss_3Cmp_fit([x',n'], k0_3);



FitPar3 = [G_Coef3;G_Sigma3];
FitPar3(1,:) = [10.^FitPar3(1,1:6),FitPar3(1,7:8)];
FitPar3(2,:) = [FitPar3(1,1:6) - 10.^(G_Coef3(1:6) - G_Sigma3(1:6)),FitPar3(2,7:8)];
FitHist3 = G_Fit3;

SSR3 = sum((DLog_hist(:,2) - G_Fit3(:,2)).^2);

R2_3 = 1 - (SSR3/SST);
FitPar3(:,9) = [SSR3;R2_3];

if R2_3 >= 0.95 && fit2Show == 0
    fit2Show = 3;
end
handles.FitPar = [FitPar1,FitPar2,FitPar3];
handles.Hist.Fit1 = FitHist1;
handles.Hist.Fit2 = FitHist2;
handles.Hist.Fit3 = FitHist3;


if fit2Show == 1
    set(handles.showFitPop,'Value',1);
    plot(G_Fit1(:,1),G_Fit1(:,2),'r');
    
    FitString{1} = ['n =', num2str(sum(handles.Data(:,2)))];
    FitString{2} = ['D1 = ', num2str(handles.FitPar(1,1),3),'+/-', num2str(handles.FitPar(2,1),3)];
    FitString{3} = ['SSR = ', num2str(handles.FitPar(1,3),3)];
    FitString{4} = ['R^2 = ', num2str(handles.FitPar(2,3),3)];
    

elseif fit2Show == 2
    set(handles.showFitPop,'Value',2);
    plot(G_Fit2(:,1),G_Fit2(:,2),'r')
    plot(G_Fit2(:,1),G_Fit2(:,3),'k--');
    plot(G_Fit2(:,1),G_Fit2(:,4),'k-.');
    
    FitString{1} = ['n =', num2str(sum(handles.Data(:,2)))];
    FitString{2} = ['D1 = ', num2str(handles.FitPar(1,4),3),'+/-', num2str(handles.FitPar(2,4),3)];
    FitString{3} = ['D2 = ', num2str(handles.FitPar(1,6),3),'+/-', num2str(handles.FitPar(2,6),3)];
    FitString{4} = ['F1 = ', num2str(handles.FitPar(1,8),3),'+/-', num2str(handles.FitPar(2,8),3)];
    FitString{5} = ['SSR = ', num2str(handles.FitPar(1,9),3)];
    FitString{6} = ['R^2 = ', num2str(handles.FitPar(2,9),3)];
else
    set(handles.showFitPop,'Value',3);
    plot(G_Fit3(:,1),G_Fit3(:,2),'r')
    plot(G_Fit3(:,1),G_Fit3(:,3),'k--');
    plot(G_Fit3(:,1),G_Fit3(:,4),'k-.');
    plot(G_Fit3(:,1),G_Fit3(:,5),'k.');
    
    FitString{1} = ['n =', num2str(sum(handles.Data(:,2)))];
    FitString{2} = ['D1 = ', num2str(handles.FitPar(1,10),3),'+/-', num2str(handles.FitPar(2,10),3)];
    FitString{3} = ['D2 = ', num2str(handles.FitPar(1,12),3),'+/-', num2str(handles.FitPar(2,12),3)];
    FitString{4} = ['D3 = ', num2str(handles.FitPar(1,14),3),'+/-', num2str(handles.FitPar(2,14),3)];
    FitString{5} = ['F1 = ', num2str(handles.FitPar(1,16),3),'+/-', num2str(handles.FitPar(2,16),3)];
    FitString{6} = ['F2 = ', num2str(handles.FitPar(1,17),3),'+/-', num2str(handles.FitPar(2,17),3)];
    FitString{7} = ['SSR = ', num2str(handles.FitPar(1,18),3)];
    FitString{8} = ['R^2 = ', num2str(handles.FitPar(2,18),3)];
end
hold off;

set(handles.FitResultsEdit,'String',FitString);
figname = get(gcf,'Name');
figname = [figname ' ' ROIChoice{1}];
set(gcf,'Name',figname);
    




% Choose default command line output for DiffCoeffMerge
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DiffCoeffMerge wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DiffCoeffMerge_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in showFitPop.
function showFitPop_Callback(hObject, eventdata, handles)
% hObject    handle to showFitPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns showFitPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from showFitPop

fit2Show = get(hObject,'Value');

axes(handles.axes1);
bar(handles.Data(:,1),handles.Data(:,2),'EdgeColor', [0 0 0], ...
'FaceColor','none', 'BarWidth', 1, 'LineWidth', 0.75);
xlabel('LogD (\mum^2/s)')
ylabel('Frequency');
hold on;

if fit2Show == 1
    set(handles.showFitPop,'Value',1);
    plot(handles.Hist.Fit1(:,1),handles.Hist.Fit1(:,2),'r');
    
    FitString{1} = ['n =', num2str(sum(handles.Data(:,2)))];
    FitString{2} = ['D1 = ', num2str(handles.FitPar(1,1),3),'+/-', num2str(handles.FitPar(2,1),3)];
    FitString{3} = ['SSR = ', num2str(handles.FitPar(1,3),3)];
    FitString{4} = ['R^2 = ', num2str(handles.FitPar(2,3),3)];
    

elseif fit2Show == 2
    set(handles.showFitPop,'Value',2);
    plot(handles.Hist.Fit2(:,1),handles.Hist.Fit2(:,2),'r')
    plot(handles.Hist.Fit2(:,1),handles.Hist.Fit2(:,3),'k--');
    plot(handles.Hist.Fit2(:,1),handles.Hist.Fit2(:,4),'k-.');
    
    FitString{1} = ['n =', num2str(sum(handles.Data(:,2)))];
    FitString{2} = ['D1 = ', num2str(handles.FitPar(1,4),3),'+/-', num2str(handles.FitPar(2,4),3)];
    FitString{3} = ['D2 = ', num2str(handles.FitPar(1,6),3),'+/-', num2str(handles.FitPar(2,6),3)];
    FitString{4} = ['F1 = ', num2str(handles.FitPar(1,8),3),'+/-', num2str(handles.FitPar(2,8),3)];
    FitString{5} = ['SSR = ', num2str(handles.FitPar(1,9),3)];
    FitString{6} = ['R^2 = ', num2str(handles.FitPar(2,9),3)];
else
    set(handles.showFitPop,'Value',3);
    plot(handles.Hist.Fit3(:,1),handles.Hist.Fit3(:,2),'r')
    plot(handles.Hist.Fit3(:,1),handles.Hist.Fit3(:,3),'k--');
    plot(handles.Hist.Fit3(:,1),handles.Hist.Fit3(:,4),'k-.');
    plot(handles.Hist.Fit3(:,1),handles.Hist.Fit3(:,5),'k:');
    
    FitString{1} = ['n =', num2str(sum(handles.Data(:,2)))];
    FitString{2} = ['D1 = ', num2str(handles.FitPar(1,10),3),'+/-', num2str(handles.FitPar(2,10),3)];
    FitString{3} = ['D2 = ', num2str(handles.FitPar(1,12),3),'+/-', num2str(handles.FitPar(2,12),3)];
    FitString{4} = ['D3 = ', num2str(handles.FitPar(1,14),3),'+/-', num2str(handles.FitPar(2,14),3)];
    FitString{5} = ['F1 = ', num2str(handles.FitPar(1,16),3),'+/-', num2str(handles.FitPar(2,16),3)];
    FitString{6} = ['F2 = ', num2str(handles.FitPar(1,17),3),'+/-', num2str(handles.FitPar(2,17),3)];
    FitString{7} = ['SSR = ', num2str(handles.FitPar(1,18),3)];
    FitString{8} = ['R^2 = ', num2str(handles.FitPar(2,18),3)];
end
hold off;
set(handles.FitResultsEdit,'String',FitString);


% --- Executes during object creation, after setting all properties.
function showFitPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to showFitPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SaveFigureButton.
function SaveFigureButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFigureButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

frameTime = handles.frameTime;
Data = handles.Data;
binSize = handles.binSize;
FitPar = handles.FitPar;

Fits = handles.Hist;

DefaultName = handles.FileName;
SearchStr = '(.*)\.\w*';

% tlist = handles.JDH(2:end, 1);
% rlist = handles.JDH(1, 2:end)';
% rlist(1,2) =tlist(1);

scrsz = get(0,'ScreenSize');
h1 = figure('Position',[.15*scrsz(3) .15*scrsz(4) .4*scrsz(3) .5*scrsz(4)]) ;

Selection = get(handles.showFitPop, 'Value');
set(gca, 'FontSize', 10);
bar(Data(:,1),Data(:,2),'EdgeColor', [0 0 0], ...
    'FaceColor','none', 'BarWidth', 1, 'LineWidth', 0.75);
hold on;

if Selection == 1 % One component diffusion fit
    
    % plot the histogram of displacements and the fit in the current
    % axes
    
    
    plot(Fits.Fit1(:,1),Fits.Fit1(:,2),'r', 'LineWidth', 1);
    hold off;
    xlabel('LogD [\mum^2/s]', 'FontSize', 10);
    ylabel('Counts', 'FontSize', 10);
    
    FitString{1} = 'One Component Diffusion Coefficient Fit';
    FitString{2} = ['n = ' num2str(sum(Data(:,2))),'; ',...
        'D = ', num2str(FitPar(1,1),3),' \pm ', num2str(FitPar(2,1),2), ' \mum^2/s; '];
    FitString{3}= ['SSR = ' num2str(FitPar(1,3)),'; \DeltaT = ', num2str(frameTime)];
    
    
    title(FitString,'FontSize', 10);
    legend({'Track Diffusion Coefficients';'One Component Fit'}, 'FontSize', 10,'Location','NorthWest');
    
    
    DefaultName = regexprep(DefaultName, SearchStr, '$1_Diffusion_1Cmp');
    FilterSpec = {'*.fig';'*.ai';'*.bpm';'*.eps';'*.jpeg';'*.pdf';'*.tif'};
    [FileNameOut,PathNameOut] = uiputfile(FilterSpec,'Save Figure',DefaultName);
    if FileNameOut ~= 0;
        saveas(h1, [PathNameOut,FileNameOut]);
    end
    
elseif Selection == 2 % Two component diffusion fit
    
    % plot the histogram of displacements and the fit in the current
    % axes
    
    plot(Fits.Fit2(:,1),Fits.Fit2(:,2),'r', 'LineWidth', 1);
    plot(Fits.Fit2(:,1),Fits.Fit2(:,3),'k--', 'LineWidth', 1);
    plot(Fits.Fit2(:,1),Fits.Fit2(:,4),'k-.', 'LineWidth', 1);
    hold off;
    xlabel('LogD [\mum^2/s]', 'FontSize', 10);
    ylabel('Counts', 'FontSize', 10);
    
    FitString{1} = 'Two Component Diffusion Coefficient Fit';
    FitString{2} = ['n = ' num2str(sum(Data(:,2)))];
    FitString{3} = ['D_1 = ', num2str(FitPar(1,4),3),' \pm ', num2str(FitPar(2,4),2), ' \mum^2/s; ',...
        'D_2 = ', num2str(FitPar(1,6),3),' \pm ', num2str(FitPar(2,6),2), ' \mum^2/s'];
    FitString{4} = ['f_1 = ', num2str(FitPar(1,8),3),' \pm ', num2str(FitPar(2,8),2)];
    FitString{5}= ['SSR = ' num2str(FitPar(1,9)),'; \DeltaT = ', num2str(frameTime)];
    
    
    title(FitString,'FontSize', 10);
    legend({'Track Diffusion Coefficients';'Two Component Fit';'Component 1'; 'Component 2'}, 'FontSize', 10,'Location','NorthWest');
    
    
    DefaultName = regexprep(DefaultName, SearchStr, '$1_Diffusion_2Cmp');
    FilterSpec = {'*.fig';'*.ai';'*.bpm';'*.eps';'*.jpeg';'*.pdf';'*.tif'};
    [FileNameOut,PathNameOut] = uiputfile(FilterSpec,'Save Figure',DefaultName);
    if FileNameOut ~= 0;
        saveas(h1, [PathNameOut,FileNameOut]);
    end
    
elseif Selection == 3 % Three component diffusion fit
    
    % plot the histogram of displacements and the fit in the current
    % axes
    
    plot(Fits.Fit3(:,1),Fits.Fit3(:,2),'r', 'LineWidth', 1);
    plot(Fits.Fit3(:,1),Fits.Fit3(:,3),'k--', 'LineWidth', 1);
    plot(Fits.Fit3(:,1),Fits.Fit3(:,4),'k-.', 'LineWidth', 1);
    plot(Fits.Fit3(:,1),Fits.Fit3(:,5),'k:', 'LineWidth', 1);
    hold off;
    xlabel('LogD [\mum^2/s]', 'FontSize', 10);
    ylabel('Counts', 'FontSize', 10);
    
    FitString{1} = 'Three Component Diffusion Coefficient Fit';
    FitString{2} = ['n = ' num2str(sum(Data(:,2)))];
    FitString{3} = ['D_1 = ', num2str(FitPar(1,10),3),' \pm ', num2str(FitPar(2,10),2), ' \mum^2/s; ',...
        'D_2 = ', num2str(FitPar(1,12),3),' \pm ', num2str(FitPar(2,12),2), ' \mum^2/s; ',...
        'D_3 = ', num2str(FitPar(1,14),3),' \pm ', num2str(FitPar(2,14),2), ' \mum^2/s'];
    FitString{4} = ['f_1 = ', num2str(FitPar(1,16),3),' \pm ', num2str(FitPar(2,16),2), '; ',...
        'f_2 = ', num2str(FitPar(1,17),3),' \pm ', num2str(FitPar(2,17),2)];
    FitString{5}= ['SSR = ' num2str(FitPar(1,18)),'; \DeltaT = ', num2str(frameTime)];
    
    
    title(FitString,'FontSize', 10);
    legend({'Track Diffusion Coefficients';'Two Component Fit';'Component 1'; 'Component 2'; 'Component 3'}, 'FontSize', 10,'Location','NorthWest');
    
    
    DefaultName = regexprep(DefaultName, SearchStr, '$1_Diffusion_3Cmp');
    FilterSpec = {'*.fig';'*.ai';'*.bpm';'*.eps';'*.jpeg';'*.pdf';'*.tif'};
    [FileNameOut,PathNameOut] = uiputfile(FilterSpec,'Save Figure',DefaultName);
    if FileNameOut ~= 0;
        saveas(h1, [PathNameOut,FileNameOut]);
    end
    
    
    
end

% --- Executes on button press in CopyFitParams.
function CopyFitParams_Callback(hObject, eventdata, handles)
% hObject    handle to CopyFitParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

num2clip(handles.FitPar);
msgbox({'The Array containing the Fit Parameters';...
    'has been copied to the clipboard'});

% --- Executes on button press in CopyHistogram.
function CopyHistogram_Callback(hObject, eventdata, handles)
% hObject    handle to CopyHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Data = handles.Data;
Data = [Data, handles.Hist.Fit1(:,2),handles.Hist.Fit2(:,2:4),handles.Hist.Fit3(:,2:5)];

num2clip(Data);
msgbox({'The Array containing the Diffusion Coefficient histogram';...
    'and all three fits has been copied to the clipboard'});
% --- Executes on button press in DoneButton.
function DoneButton_Callback(hObject, eventdata, handles)
% hObject    handle to DoneButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close();

function FitResultsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to FitResultsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FitResultsEdit as text
%        str2double(get(hObject,'String')) returns contents of FitResultsEdit as a double


% --- Executes during object creation, after setting all properties.
function FitResultsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FitResultsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
