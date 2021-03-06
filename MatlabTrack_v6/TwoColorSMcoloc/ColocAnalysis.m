function varargout = ColocAnalysis(varargin)
% COLOCANALYSIS MATLAB code for ColocAnalysis.fig
%      COLOCANALYSIS, by itself, creates a new COLOCANALYSIS or raises the existing
%      singleton*.
%
%      H = COLOCANALYSIS returns the handle to a new COLOCANALYSIS or the handle to
%      the existing singleton*.
%
%      COLOCANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COLOCANALYSIS.M with the given input arguments.
%
%      COLOCANALYSIS('Property','Value',...) creates a new COLOCANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ColocAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ColocAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ColocAnalysis

% Last Modified by GUIDE v2.5 09-Dec-2016 16:06:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ColocAnalysis_OpeningFcn, ...
    'gui_OutputFcn',  @ColocAnalysis_OutputFcn, ...
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


% --- Executes just before ColocAnalysis is made visible.
function ColocAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ColocAnalysis (see VARARGIN)

% Choose default command line output for ColocAnalysis
handles.output = hObject;
AnalysisParam = varargin{1};
handles.Particles1 = varargin{2};
handles.Particles2 = varargin{3};
handles.Tracks1 = varargin{4};
handles.Tracks2 = varargin{5};
handles.Parameters.ThreshL = str2double(AnalysisParam{1});
% handles.Parameters.ThreshH = str2num(AnalysisS{2});
handles.Parameters.minBoundFrames = str2double(AnalysisParam{2});
handles.Parameters.bin = str2double(AnalysisParam{3});
handles.Parameters.frameTime = str2double(AnalysisParam{4});

%update the displayed analysis parameters
set(handles.BoundThreshLedit,'String',AnalysisParam{1});
set(handles.minBoundEdit,'String',AnalysisParam{2});
set(handles.HistBinEdit,'String',AnalysisParam{3});
set(handles.frameTimeEdit,'String',AnalysisParam{4});


nTimes = max(max(handles.Particles1(:,6)),max(handles.Particles2(:,6)));
TimePoints = (0:nTimes-1)';
CumParticles1 = zeros(nTimes,2);
CumParticles2 = zeros(nTimes,2);

for i = 1:nTimes
    tmp1 = handles.Particles1(handles.Particles1(:,6) == i,:);
    tmp2 = handles.Particles2(handles.Particles2(:,6) == i,:);
    CumParticles1(i,2) = size(tmp1,1);
    CumParticles2(i,2) = size(tmp2,1);
    
end

CumParticles1(:,1) = TimePoints*handles.Parameters.frameTime;
CumParticles2(:,1) = TimePoints*handles.Parameters.frameTime;

[closeParticles, ~] = calculateColocParticles(handles.Particles1,...
    handles.Particles2,handles.Parameters.ThreshL);

PartOrTracks = questdlg('Do you want to consider all Particles or only those that were Tracked?','Choose Particles','All','Tracked','All');

% for i = 1:max(closeParticles(:,5))
%     curTrack = closeParticles(closeParticles(:,5) == i,:);
%     totalCloseParticles = totalCloseParticles + size(curTrack,1);
% end

[boundParticles, allboundParticles] = calculateColocTracks(handles.Tracks1,...
    handles.Tracks2,handles.Parameters.ThreshL,handles.Parameters.minBoundFrames);
if strcmp(PartOrTracks,'All')
    totalCloseParticles = size(closeParticles,1);
else
    totalCloseParticles = size(allboundParticles,1);
end


if isempty(boundParticles)
    msgbox('No Colocalized particles exist in the current image using the specified parameters','No Colocalization');
    
    close('ColocAnalysis');
    handles.output = [];
    return;
end

if isempty(boundParticles)
    set(handles.k1_1,'String','N/A');
    set(handles.SSR1,'String','N/A');
    set(handles.k1_2,'String','N/A');
    set(handles.k2_2,'String','N/A');
    set(handles.F1_2,'String','N/A');
    set(handles.SSR2,'String','N/A');
    msgbox('No data matches the criteria','No Colocalized Particles Found');
    
else

    %Calculate the duration of colocalization events
    handles.Colocalization.Duration = zeros(max(boundParticles(:,5)),1);
    for i = 1:max(boundParticles(:,5))
        CurPart = boundParticles(boundParticles(:,5) == i,:);
        handles.Colocalization.Duration(i,:) = (CurPart(end,1) - CurPart(1,1))+1;
    end


    binC = handles.Parameters.minBoundFrames:handles.Parameters.bin:max(handles.Colocalization.Duration);
    binC = (binC - 1)*handles.Parameters.frameTime;
    handles.Colocalization.Duration = (handles.Colocalization.Duration - 1)*handles.Parameters.frameTime;
    if length(binC) == 1
        x = binC;
        n = size(handles.Colocalization,1);
        N = n;
        resTimes = N;
    else
        [n,x] = hist(handles.Colocalization.Duration, binC);
        n2 = flip(n);
        
        N2 = cumsum(n2);%./sum(n);
        N = flip(N2);
        resTimes = N(1:end-1) - N(2:end);
        resTimes = [resTimes, N(end)];
    end



    answer = questdlg('Do you want to correct the histogram for photobleaching?');

    if strcmp(answer,'Yes')
        disp('____________________________________')
        disp('Estimating bleaching characteristics for Channel 1')
        disp('____________________________________')

        [BleachRates1, Dummy, CumParticles1(:,3:6)] =...
            ExpDecay_2Cmp_fit(CumParticles1, [1 0.1]);
        disp(['Bleach Rate 1: ', num2str(BleachRates1(1), 3), ' s^-1'])
        disp(['Bleach Rate 2: ', num2str(BleachRates1(2), 3), ' s^-1'])
        disp(['Fraction 1: ', num2str(BleachRates1(3), 3)])
        disp('____________________________________')
        disp(' ')

        disp('____________________________________')
        disp('Estimating bleaching characteristics for Channel 2')
        disp('____________________________________')

        [BleachRates2, Dummy, CumParticles2(:,3:6)] =...
            ExpDecay_2Cmp_fit(CumParticles2, [1 0.1]);
        disp(['Bleach Rate 1: ', num2str(BleachRates2(1), 3), ' s^-1'])
        disp(['Bleach Rate 2: ', num2str(BleachRates2(2), 3), ' s^-1'])
        disp(['Fraction 1: ', num2str(BleachRates2(3), 3)])
        disp('____________________________________')
        disp(' ')
        % Recalculate the survival probability corrected for photobleaching
        N = N./ ...
            (BleachRates1(3)*exp(-BleachRates1(1).* x) + ...
            (1-BleachRates1(3))*exp(-BleachRates1(2).* x));
        N = N./...
            (BleachRates2(3)*exp(-BleachRates2(1).* x) + ...
            (1-BleachRates2(3))*exp(-BleachRates2(2).* x));
%         resTimes = resTimes./...
%             (BleachRates2(3)*exp(-BleachRates2(1).* x) + ...
%             (1-BleachRates2(3))*exp(-BleachRates2(2).* x));
        

        set(handles.PhotoBleachCheck,'Value',1);
    else
        set(handles.PhotoBleachCheck,'Value',0);

    end

%     axes(handles.PlotAxis);
%     hold off
%     bar(x,N./N(1));

    

    Data = [x', N'];
    resTimes = [x', resTimes'];
    
    normQuest = questdlg('Do you want to Normalize to the Fraction of Proximal Molecules, or to the Fraction of All Molecules?','Normalization','Proximal','All','No Normalization','All');
    if strcmp(normQuest,'All')
        % Then calculate the total number of bound spots
        TotalBoundMolecules = 2*sum(resTimes(:,2).*((resTimes(:,1)/handles.Parameters.frameTime)+1));
        
        % and the total number of spots
        if strcmp(PartOrTracks,'All')
            TotalMolecules = sum(CumParticles1(:,2))+sum(CumParticles2(:,2));
        else
            TotalMolecules = size(handles.Tracks1,1) + size(handles.Tracks2,1);
        end
        
        
        % Finally divide the two.
        PartialBoundFraction = TotalBoundMolecules/TotalMolecules;
        
        BFerror = (sqrt(TotalBoundMolecules)/TotalBoundMolecules + ...
            sqrt(TotalMolecules)/TotalMolecules)*PartialBoundFraction;
        
        % Normalize the Survival probability for the partial bound fraction
        Data(:,2) = Data(:,2)/Data(1,2)...
            .*PartialBoundFraction;
        set(handles.AllNorm,'Value',1);
        set(handles.NoNorm,'Value',0);
        set(handles.CloseNorm,'Value',0);
    elseif strcmp(normQuest,'Proximal')
        TotalBoundMolecules = 2*sum(resTimes(:,2).*((resTimes(:,1)/handles.Parameters.frameTime)+1));
        
        TotalMolecules = 2*totalCloseParticles;
        PartialBoundFraction = TotalBoundMolecules/TotalMolecules;
        BFerror = (sqrt(TotalBoundMolecules)/TotalBoundMolecules + ...
            sqrt(TotalMolecules)/TotalMolecules)*PartialBoundFraction;
        Data(:,2) = Data(:,2)/Data(1,2)...
            .*PartialBoundFraction;
        set(handles.AllNorm,'Value',0);
        set(handles.NoNorm,'Value',0);
        set(handles.CloseNorm,'Value',1);
    else
        Data(:,2) = Data(:,2)/Data(1,2);
        set(handles.AllNorm,'Value',0);
        set(handles.NoNorm,'Value',1);
        set(handles.CloseNorm,'Value',0);
    end
    handles.TotalMolecules = TotalMolecules;
    handles.TotalBoundMolecules = TotalBoundMolecules;
    % 
    % % NOTE: THIS IS ONLY THE PARTIAL BOUND FRACTION
    % % BECAUSE WE ARE LOOKING ONLY AT PARTICLES BOUND
    % % FOR MORE THAN Nmin. We HAVE TO FIT THE HISTOGRAM OF RESIDENCE
    % % TIMES AND EXTRAPOLATE THE TRUE VALUE OF THE BOUND FRACTION
    axes(handles.PlotAxis);
    hold off
    bar(Data(:,1),Data(:,2));
    
    xlabel('Time (seconds)');
    ylabel('Fraction of counts');
    if size(Data,1) < 2
        set(handles.k1_1,'String','N/A');
        set(handles.SSR1,'String','N/A');
        set(handles.k1_2,'String','N/A');
        set(handles.k2_2,'String','N/A');
        set(handles.F1_2,'String','N/A');
        set(handles.SSR2,'String','N/A');
        
        handles.Colocalization.Hist = Data;
        handles.Colocalization.OneCompFit = [];
        handles.Colocalization.TwoCompFit = [];
        handles.FitPar = [];
        msgbox('Not enough points in the histogram to perform fitting','1-component fit');

    elseif size(Data,1) < 4
        % Fit exponentials to the Survival probability
        [fitpar,espSigma, fit]= ExpDecay_fit(Data, 1/(Data(1)*10));
        hold on;
        plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
        legend(['Colocalization Times (C_{eq} = ', num2str(Data(1,2),3),')'], ['Single Exponential Fit (C_{eq} = ',num2str(fitpar(2),3),')']);
        set(handles.k1_1,'String',[num2str(fitpar(1),3),'+/-',num2str(espSigma(1),2)]);
        SSR1 = sum((fit(:,2)-Data(:,2)).^2);
        set(handles.SSR1,'String',num2str(SSR1,4));
        handles.FitPar = [fitpar; espSigma];
        set(handles.k1_2,'String','N/A');
        set(handles.k2_2,'String','N/A');
        set(handles.F1_2,'String','N/A');
        set(handles.SSR2,'String','N/A');
        
        handles.Colocalization.Hist = Data;
        handles.Colocalization.OneCompFit = fit;
        handles.Colocalization.TwoCompFit = [];
        handles.FitPar = [fitpar, zeros(1,4); espSigma, zeros(1,4)];

        msgbox('Not enough points in histogram to fit with 2-component model','2-component fit');

    else
        % Fit exponentials to the Survival probability
        [fitpar,espSigma, fit]= ExpDecay_fit(Data, 1/(Data(1)*10));

        [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit(Data, [1/(Data(1)*10) 1/(Data(1)*100)]);
        Data = [Data,fit(:,2),fit2(:,2:4)];
        hold on;

        plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
        plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
        plot(fit2(:,1),fit2(:,3),'c', 'LineWidth', 1);
        plot(fit2(:,1),fit2(:,4),'k', 'LineWidth', 1);

        legend(['Colocalization Times (C_{eq} = ', num2str(Data(1,2),3),')'], ['Single Exponential Fit (C_{eq} = ',num2str(fitpar(2),3),')'], ['Double Exponential Fit (C_{eq} = ',num2str(fitpar2(4),3),')'],'Component 1','Component 2');

        %Copy the fitting parameters to the appropriate text boxes
        set(handles.k1_1,'String',[num2str(fitpar(1),3),' +/- ',num2str(espSigma(1),2)]);

        set(handles.k1_2,'String',[num2str(fitpar2(1),3),' +/- ',num2str(espSigma2(1),2)]);
        set(handles.k2_2,'String',[num2str(fitpar2(2),3),' +/- ',num2str(espSigma2(2),2)]);
        set(handles.F1_2,'String',[num2str(fitpar2(3),2),' +/- ',num2str(espSigma2(3),2)]);

        SSR1 = sum((fit(:,2)-Data(:,2)).^2);

        SSR2 = sum((fit2(:,2)-Data(:,2)).^2);
        
        [~, pvalRes,FstatRes] = FtestModelCompare(Data(:,2),fit(:,2),fit2(:,2),2,4);
        set(handles.pValEdit,'String',num2str(pvalRes));
        set(handles.FstatEdit,'String',num2str(FstatRes));
        
        
        set(handles.SSR1,'String',num2str(SSR1,4));
        set(handles.SSR2,'String',num2str(SSR2,4));
        handles.Colocalization.Hist = Data;
        handles.Colocalization.OneCompFit = fit;
        handles.Colocalization.TwoCompFit = fit2;
        handles.FitPar = [fitpar, fitpar2, pvalRes; espSigma, espSigma2,FstatRes];
    end
    set(handles.BoundMolText,'String',num2str(handles.TotalBoundMolecules));
    set(handles.TotalMolText,'String',num2str(handles.TotalMolecules));
    
    
    set(handles.Tracks1Text,'String',num2str(max(handles.Tracks1(:,4))));
    set(handles.Tracks2Text,'String',num2str(max(handles.Tracks2(:,4))));
   
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ColocAnalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ColocAnalysis_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = [];


% --- Executes on button press in Graph2Image.
function Graph2Image_Callback(hObject, eventdata, handles)
% hObject    handle to Graph2Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filterSpec = {'*.tif'; '*.jpg'};

defName = 'Colocalization Times';

[fname, pname,filterIndex] = uiputfile(filterSpec,'Save Graph as Image',defName);
if fname ~= 0
    scrsz = get(0,'ScreenSize');
    h1 = figure('Position',[.15*scrsz(3) .15*scrsz(4) .5*scrsz(3) .5*scrsz(4)]) ;
    Data = handles.Colocalization.Hist;
%     bar(Data(:,1),Data(:,2)./Data(1,2));
    bar(Data(:,1),Data(:,2));
    xlabel('Time (seconds)');
    ylabel('Fraction of counts');
    hold all
    fit = handles.Colocalization.OneCompFit;
    fit2 = handles.Colocalization.TwoCompFit;
%     plot(fit(:,1),fit(:,2)./fit(1,2),'g', 'LineWidth', 1);
%     plot(fit2(:,1),fit2(:,2)./fit2(1,2),'r', 'LineWidth', 1);
%     plot(fit2(:,1),fit2(:,3)./fit2(1,3),'c', 'LineWidth', 1);
%     plot(fit2(:,1),fit2(:,4)./fit2(1,4),'k', 'LineWidth', 1);
    plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,3),'c', 'LineWidth', 1);
    plot(fit2(:,1),fit2(:,4),'k', 'LineWidth', 1);
    legend('Colocalization Times', 'Single exponential fit', 'Double Exponential Fit','Component 1','Component 2');
    
    if filterIndex == 1
        outspec = '-dtiff';
    else
        outspec = '-djpeg';
    end
    print(h1,[pname fname],outspec,'-r300');
    close(h1);

end

guidata(hObject,handles);

% --- Executes on button press in HistData2Clip.
function HistData2Clip_Callback(hObject, eventdata, handles)
% hObject    handle to HistData2Clip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

num2clip(handles.Colocalization.Hist);
msgbox('Histogram Data copied to clipboard');


% --- Executes on button press in Param2Clip.
function Param2Clip_Callback(hObject, eventdata, handles)
% hObject    handle to Param2Clip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

num2clip(handles.FitPar);

msgbox({'Fit Parameters copied to clipboard';...
    'Format: {k1_1}   {A0_1}    {k1_2}    {k2_2}    {F1_2}    {A0_2} {F p-val}' ;...
    '      {epsk1_1} {epsA0_1} {epsk1_2} {epsk2_2} {epsF1_2} {epsA0_2} {F stat}'});

function BoundThreshLedit_Callback(hObject, eventdata, handles)
% hObject    handle to BoundThreshLedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BoundThreshLedit as text
%        str2double(get(hObject,'String')) returns contents of BoundThreshLedit as a double


% --- Executes during object creation, after setting all properties.
function BoundThreshLedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BoundThreshLedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minBoundEdit_Callback(hObject, eventdata, handles)
% hObject    handle to minBoundEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minBoundEdit as text
%        str2double(get(hObject,'String')) returns contents of minBoundEdit as a double


% --- Executes during object creation, after setting all properties.
function minBoundEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minBoundEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function HistBinEdit_Callback(hObject, eventdata, handles)
% hObject    handle to HistBinEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HistBinEdit as text
%        str2double(get(hObject,'String')) returns contents of HistBinEdit as a double


% --- Executes during object creation, after setting all properties.
function HistBinEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HistBinEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PhotoBleachCheck.
function PhotoBleachCheck_Callback(hObject, eventdata, handles)
% hObject    handle to PhotoBleachCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PhotoBleachCheck


% --- Executes on button press in AnalyzePush.
function AnalyzePush_Callback(hObject, eventdata, handles)
% hObject    handle to AnalyzePush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Parameters.ThreshL = str2double(get(handles.BoundThreshLedit,'String'));

handles.Parameters.minBoundFrames = str2double(get(handles.minBoundEdit,'String'));
handles.Parameters.bin = str2double(get(handles.HistBinEdit,'String'));
handles.Parameters.frameTime = str2double(get(handles.frameTimeEdit,'String'));

nTimes = max(max(handles.Particles1(:,6)),max(handles.Particles2(:,6)));
TimePoints = (0:nTimes-1)';
CumParticles1 = zeros(nTimes,2);
CumParticles2 = zeros(nTimes,2);


for i = 1:nTimes
    tmp1 = handles.Particles1(handles.Particles1(:,6) == i,:);
    tmp2 = handles.Particles2(handles.Particles2(:,6) == i,:);
    CumParticles1(i,2) = size(tmp1,1);
    CumParticles2(i,2) = size(tmp2,1);
    
end

CumParticles1(:,1) = TimePoints*handles.Parameters.frameTime;
CumParticles2(:,1) = TimePoints*handles.Parameters.frameTime;



[closeParticles, ~] = calculateColocParticles(handles.Particles1,...
    handles.Particles2,handles.Parameters.ThreshL);
% totalCloseParticles = size(closeParticles,1);
% for i = 1:max(closeParticles(:,5))
%     curTrack = closeParticles(closeParticles(:,5) == i,:);
%     totalCloseParticles = totalCloseParticles + size(curTrack,1);
% end



[boundParticles, allboundParticles] = calculateColocTracks(handles.Tracks1,...
    handles.Tracks2,handles.Parameters.ThreshL,handles.Parameters.minBoundFrames);
PartOrTracks = questdlg('Do you want to consider all Particles or only those that were Tracked?','Choose Particles','All','Tracked','All');

if strcmp(PartOrTracks,'All')
    totalCloseParticles = size(closeParticles,1);
else
    totalCloseParticles = size(allboundParticles,1);
end
if isempty(boundParticles)
    set(handles.k1_1,'String','N/A');
    set(handles.SSR1,'String','N/A');
    set(handles.k1_2,'String','N/A');
    set(handles.k2_2,'String','N/A');
    set(handles.F1_2,'String','N/A');
    set(handles.SSR2,'String','N/A');
    msgbox('No data matches the criteria','No Colocalized Particles Found');
    
else

    %Calculate the duration of colocalization events
    handles.Colocalization.Duration = zeros(max(boundParticles(:,5)),1);
    for i = 1:max(boundParticles(:,5))
        CurPart = boundParticles(boundParticles(:,5) == i,:);
        handles.Colocalization.Duration(i,:) = (CurPart(end,1) - CurPart(1,1))+1;
    end


    binC = handles.Parameters.minBoundFrames:handles.Parameters.bin:max(handles.Colocalization.Duration);
    binC = (binC - 1)*handles.Parameters.frameTime;
    
    handles.Colocalization.Duration = (handles.Colocalization.Duration - 1)*handles.Parameters.frameTime;
    if length(binC) == 1
        x = binC;
        n = size(handles.Colocalization,1);
        N = n;
        resTimes = N;
    else
        [n,x] = hist(handles.Colocalization.Duration, binC);
        n2 = flip(n);
        N2 = cumsum(n2);%./sum(n);
        N = flip(N2);
        resTimes = N(1:end-1) - N(2:end);
        resTimes = [resTimes, N(end)];
        
    end
    



    if get(handles.PhotoBleachCheck,'Value')
        disp('____________________________________')
        disp('Estimating bleaching characteristics for Channel 1')
        disp('____________________________________')

        [BleachRates1, Dummy, CumParticles1(:,3:6)] =...
            ExpDecay_2Cmp_fit(CumParticles1, [1 0.1]);
        disp(['Bleach Rate 1: ', num2str(BleachRates1(1), 3), ' s^-1'])
        disp(['Bleach Rate 2: ', num2str(BleachRates1(2), 3), ' s^-1'])
        disp(['Fraction 1: ', num2str(BleachRates1(3), 3)])
        disp('____________________________________')
        disp(' ')

        disp('____________________________________')
        disp('Estimating bleaching characteristics for Channel 2')
        disp('____________________________________')

        [BleachRates2, Dummy, CumParticles2(:,3:6)] =...
            ExpDecay_2Cmp_fit(CumParticles2, [1 0.1]);
        disp(['Bleach Rate 1: ', num2str(BleachRates2(1), 3), ' s^-1'])
        disp(['Bleach Rate 2: ', num2str(BleachRates2(2), 3), ' s^-1'])
        disp(['Fraction 1: ', num2str(BleachRates2(3), 3)])
        disp('____________________________________')
        disp(' ')
        % Recalculate the survival probability corrected for photobleaching
        N = N./ ...
            (BleachRates1(3)*exp(-BleachRates1(1).* x) + ...
            (1-BleachRates1(3))*exp(-BleachRates1(2).* x));
        N = N./...
            (BleachRates2(3)*exp(-BleachRates2(1).* x) + ...
            (1-BleachRates2(3))*exp(-BleachRates2(2).* x));
%         resTimes = resTimes./...
%             (BleachRates2(3)*exp(-BleachRates2(1).* x) + ...
%             (1-BleachRates2(3))*exp(-BleachRates2(2).* x));

        set(handles.PhotoBleachCheck,'Value',1);
    else
        set(handles.PhotoBleachCheck,'Value',0);

    end

    
    Data = [x', N']; %CDF
    resTimes = [x', resTimes'];
    
    
    
    if get(handles.AllNorm,'Value')
        % Then calculate the total number of bound spots
        TotalBoundMolecules = 2*sum(resTimes(:,2).*((resTimes(:,1)/handles.Parameters.frameTime)+1));
        
        % and the total number of spots
        if strcmp(PartOrTracks,'All')
            TotalMolecules = sum(CumParticles1(:,2))+sum(CumParticles2(:,2));
        else
            TotalMolecules = size(handles.Tracks1,1) + size(handles.Tracks2,1);
        end
%         TotalMolecules = size(handles.Tracks1,1) + size(handles.Tracks2,1);
        % Finally divide the two.
        PartialBoundFraction = TotalBoundMolecules/TotalMolecules;
        
        BFerror = (sqrt(TotalBoundMolecules)/TotalBoundMolecules + ...
            sqrt(TotalMolecules)/TotalMolecules)*PartialBoundFraction;
        
        % Normalize the Survival probability for the partial bound fraction
        Data(:,2) = Data(:,2)/Data(1,2)...
            .*PartialBoundFraction;
        set(handles.AllNorm,'Value',1);
        set(handles.NoNorm,'Value',0);
        set(handles.CloseNorm,'Value',0);
    elseif get(handles.CloseNorm,'Value')
        TotalBoundMolecules = 2*sum(resTimes(:,2).*((resTimes(:,1)/handles.Parameters.frameTime)+1));
        
        TotalMolecules = 2*totalCloseParticles;%/handles.Parameters.frameTime;
        PartialBoundFraction = TotalBoundMolecules/TotalMolecules;
        BFerror = (sqrt(TotalBoundMolecules)/TotalBoundMolecules + ...
            sqrt(TotalMolecules)/TotalMolecules)*PartialBoundFraction;
        Data(:,2) = Data(:,2)/Data(1,2)...
            .*PartialBoundFraction;
        set(handles.AllNorm,'Value',0);
        set(handles.NoNorm,'Value',0);
        set(handles.CloseNorm,'Value',1);
    else
        Data(:,2) = Data(:,2)/Data(1,2);
        set(handles.AllNorm,'Value',0);
        set(handles.NoNorm,'Value',1);
        set(handles.CloseNorm,'Value',0);
    end
    handles.TotalMolecules = TotalMolecules;
    handles.TotalBoundMolecules = TotalBoundMolecules;
    
    % NOTE: THIS IS ONLY THE PARTIAL BOUND FRACTION
    % BECAUSE WE ARE LOOKING ONLY AT PARTICLES BOUND
    % FOR MORE THAN Nmin. We HAVE TO FIT THE HISTOGRAM OF RESIDENCE
    % TIMES AND EXTRAPOLATE THE TRUE VALUE OF THE BOUND FRACTION
    axes(handles.PlotAxis);
    hold off
    bar(Data(:,1),Data(:,2));
    
    xlabel('Time (seconds)');
    ylabel('Fraction of counts');

    if size(Data,1) < 2
        set(handles.k1_1,'String','N/A');
        set(handles.SSR1,'String','N/A');
        set(handles.k1_2,'String','N/A');
        set(handles.k2_2,'String','N/A');
        set(handles.F1_2,'String','N/A');
        set(handles.SSR2,'String','N/A');
        msgbox('Not enough points in the histogram to perform fitting','1-component fit');

    elseif size(Data,1) < 4
        % Fit exponentials to the Survival probability
        [fitpar,espSigma, fit]= ExpDecay_fit(Data, min(0.5,1/(Data(1)*10)));
        hold on;
        plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
        legend(['Colocalization Times (C_{eq} = ', num2str(Data(1,2),3),')'], ['Single Exponential Fit (C_{eq} = ',num2str(fitpar(2),3),')']);
        set(handles.k1_1,'String',[num2str(fitpar(1),3),' +/- ',num2str(espSigma(1),2)]);
        SSR1 = sum((fit(:,2)-Data(:,2)).^2);
        set(handles.SSR1,'String',num2str(SSR1,4));
        handles.FitPar = [fitpar; espSigma];
        set(handles.k1_2,'String','N/A');
        set(handles.k2_2,'String','N/A');
        set(handles.F1_2,'String','N/A');
        set(handles.SSR2,'String','N/A');
        
        handles.Colocalization.Hist = Data;
        handles.Colocalization.OneCompFit = fit;
        handles.Colocalization.TwoCompFit = [];
        handles.FitPar = [fitpar, zeros(1,4); espSigma, zeros(1,4)];

        msgbox('Not enough points in histogram to fit with 2-component model','2-component fit');

    else
        % Fit exponentials to the Survival probability
        [fitpar,espSigma, fit]= ExpDecay_fit(Data, min(0.5,1/(Data(1)*10)));

        [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit(Data, [min(0.5,1/(Data(1)*10)) min(0.05,1/(Data(1)*100))]);
        Data = [Data,fit(:,2),fit2(:,2:4)];
        
        hold on;

%         plot(fit(:,1),fit(:,2)./fit(1,2),'g', 'LineWidth', 1);
        plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
        
%         plot(fit2(:,1),fit2(:,2)./fit2(1,2),'r', 'LineWidth', 1);
        plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
         
%         plot(fit2(:,1),fit2(:,3)./fit2(1,3),'c', 'LineWidth', 1);
        plot(fit2(:,1),fit2(:,3),'c', 'LineWidth', 1);
        
%         plot(fit2(:,1),fit2(:,4)./fit2(1,4),'k', 'LineWidth', 1);
        plot(fit2(:,1),fit2(:,4),'k', 'LineWidth', 1);

        legend(['Colocalization Times (C_{eq} = ', num2str(Data(1,2),3),')'], ['Single Exponential Fit (C_{eq} = ',num2str(fitpar(2),3),')'], ['Double Exponential Fit (C_{eq} = ',num2str(fitpar2(4),3),')'],'Component 1','Component 2');
        
        [~, pvalRes,FstatRes] = FtestModelCompare(Data(:,2),fit(:,2),fit2(:,2),2,4);
        set(handles.pValEdit,'String',num2str(pvalRes));
        set(handles.FstatEdit,'String',num2str(FstatRes));
        
        %Copy the fitting parameters to the appropriate text boxes
        
        
        set(handles.k1_1,'String',[num2str(fitpar(1),3),' +/- ',num2str(espSigma(1),2)]);
        set(handles.k1_2,'String',[num2str(fitpar2(1),3),' +/- ',num2str(espSigma2(1),2)]);
        set(handles.k2_2,'String',[num2str(fitpar2(2),3),' +/- ',num2str(espSigma2(2),2)]);
        set(handles.F1_2,'String',[num2str(fitpar2(3),2),' +/- ',num2str(espSigma2(3),2)]);

        SSR1 = sum((fit(:,2)-Data(:,2)).^2);

        SSR2 = sum((fit2(:,2)-Data(:,2)).^2);
        set(handles.SSR1,'String',num2str(SSR1,4));
        set(handles.SSR2,'String',num2str(SSR2,4));
        handles.Colocalization.Hist = Data;
        handles.Colocalization.OneCompFit = fit;
        handles.Colocalization.TwoCompFit = fit2;
        handles.FitPar = [fitpar, fitpar2, pvalRes; espSigma, espSigma2, FstatRes];
    end
    set(handles.BoundMolText,'String',num2str(handles.TotalBoundMolecules));
    set(handles.TotalMolText,'String',num2str(handles.TotalMolecules));
   
   
    set(handles.Tracks1Text,'String',num2str(max(handles.Tracks1(:,4))));
    set(handles.Tracks2Text,'String',num2str(max(handles.Tracks2(:,4))));
    
end
% Update handles structure
guidata(hObject, handles);
% 


function frameTimeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to frameTimeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frameTimeEdit as text
%        str2double(get(hObject,'String')) returns contents of frameTimeEdit as a double


% --- Executes during object creation, after setting all properties.
function frameTimeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameTimeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NoNorm.
function NoNorm_Callback(hObject, eventdata, handles)
% hObject    handle to NoNorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of NoNorm
if get(hObject,'Value')
    set(handles.CloseNorm,'Value',0);
    set(handles.AllNorm,'Value',0);
else
    set(handles.CloseNorm,'Value',1);
end

% --- Executes on button press in CloseNorm.
function CloseNorm_Callback(hObject, eventdata, handles)
% hObject    handle to CloseNorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CloseNorm
if get(hObject,'Value')
    set(handles.NoNorm,'Value',0);
    set(handles.AllNorm,'Value',0);
else
    set(handles.NoNorm,'Value',1);
end

% --- Executes on button press in AllNorm.
function AllNorm_Callback(hObject, eventdata, handles)
% hObject    handle to AllNorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AllNorm
if get(hObject,'Value')
    set(handles.CloseNorm,'Value',0);
    set(handles.NoNorm,'Value',0);
else
    set(handles.NoNorm,'Value',1);
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SaveMenu_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fname,pname] = uiputfile('*.mat','Save Colocalization Analysis','ColocAnalysis.mat');

if fname ~=0
    Results.Hist = handles.Colocalization.Hist;
    Results.FitPar = handles.FitPar;
    Results.Parameters = handles.Parameters;
    Results.TotalMolecules = handles.TotalMolecules; 
    Results.TotalBoundMolecules = handles.TotalBoundMolecules;
    Results.Tracks1 = handles.Tracks1;
    Results.Tracks2 = handles.Tracks2;
    save(fullfile(pname,fname),'Results');
    
end


% --------------------------------------------------------------------
function LoadMenu_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,pname] = uigetfile('*.mat','Load Colocalization Analysis');

if fname ~= 0
    Temp = load(fullfile(pname,fname));
    %update the appropriate variables in handles
    handles.Colocalization.Hist = Temp.Results.Hist;
    handles.FitPar = Temp.Results.FitPar;
    handles.Parameters = Temp.Results.Parameters;
    handles.TotalMolecules = Temp.Results.TotalMolecules; 
    handles.TotalBoundMolecules = Temp.Results.TotalBoundMolecules;
    if isfield(Temp.Results,'Tracks1')
        handles.Tracks1 = Temp.Results.Tracks1;
        handles.Tracks2 = Temp.Results.Tracks2;
    end
    
    %plot the data
    axes(handles.PlotAxis);
    hold off
    bar(handles.Colocalization.Hist(:,1),handles.Colocalization.Hist(:,2));
    
    xlabel('Time (seconds)');
    ylabel('Fraction of counts');
    hold on;

      
    plot(handles.Colocalization.Hist(:,1),handles.Colocalization.Hist(:,3),'g', 'LineWidth', 1);
    
    
    plot(handles.Colocalization.Hist(:,1),handles.Colocalization.Hist(:,4),'r', 'LineWidth', 1);
    
    
    plot(handles.Colocalization.Hist(:,1),handles.Colocalization.Hist(:,5),'c', 'LineWidth', 1);
    
    
    plot(handles.Colocalization.Hist(:,1),handles.Colocalization.Hist(:,6),'k', 'LineWidth', 1);
    
    legend(['Colocalization Times (C_{eq} = ', num2str(handles.Colocalization.Hist(1,2),3),')'], ...
        ['Single Exponential Fit (C_{eq} = ',num2str(handles.FitPar(1,2),3),')'], ...
        ['Double Exponential Fit (C_{eq} = ',num2str(handles.FitPar(1,6),3),')'],...
        'Component 1','Component 2');
   %update the fitting parameter displays  
   set(handles.k1_1,'String',[num2str(handles.FitPar(1,1),3),' +/- ',num2str(handles.FitPar(2,1),2)]);
   set(handles.k1_2,'String',[num2str(handles.FitPar(1,3),3),' +/- ',num2str(handles.FitPar(2,3),2)]);
   set(handles.k2_2,'String',[num2str(handles.FitPar(1,4),3),' +/- ',num2str(handles.FitPar(2,4),2)]);
   set(handles.F1_2,'String',[num2str(handles.FitPar(1,5),2),' +/- ',num2str(handles.FitPar(2,5),2)]);
   
   SSR1 = sum((handles.Colocalization.Hist(:,3)-handles.Colocalization.Hist(:,2)).^2);
   
   SSR2 = sum((handles.Colocalization.Hist(:,4)-handles.Colocalization.Hist(:,2)).^2);
   set(handles.SSR1,'String',num2str(SSR1,4));
   set(handles.SSR2,'String',num2str(SSR2,4));
   
   %update the analysis parameters
   set(handles.BoundThreshLedit,'String',num2str(handles.Parameters.ThreshL));
   set(handles.minBoundEdit,'String',num2str(handles.Parameters.minBoundFrames));
   set(handles.HistBinEdit,'String',num2str(handles.Parameters.bin));
   set(handles.frameTimeEdit,'String',num2str(handles.Parameters.frameTime));
   
   set(handles.BoundMolText,'String',num2str(handles.TotalBoundMolecules));
   set(handles.TotalMolText,'String',num2str(handles.TotalMolecules));
   
   if isfield(Temp.Results,'Tracks1')
       set(handles.Tracks1Text,'String',num2str(max(handles.Tracks1(:,4))));
       set(handles.Tracks2Text,'String',num2str(max(handles.Tracks2(:,4))));
   end
   
   
   
    guidata(hObject,handles);
end
    
    
