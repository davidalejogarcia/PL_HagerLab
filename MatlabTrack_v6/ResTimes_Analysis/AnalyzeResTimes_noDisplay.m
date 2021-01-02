function [Hist, FitPar] = AnalyzeResTimes_noDisplay(fileName,Parameters,BleachCurveFromInt)

Temp = load(fileName,'Results');
if strcmp(fileName,'13_ JF549 H2B-Halo_t350_300tracks substack_7.mat');
    nao = 10;
end
if isfield(Temp.Results, 'PreAnalysis') && ...
        isfield(Temp.Results.PreAnalysis, 'Tracks_um');
    Tracks = Temp.Results.PreAnalysis.Tracks_um;
    NParticles = Temp.Results.PreAnalysis.NParticles;
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
                    Tracks = Tracks(Tracks(:,5) == ROIidx,:);
                    NParticles = [NParticles(:,1) NParticles(:,ROIidx+1)];
                    ROIChoice{end+1,1} = Temp.Results.Process.ROIlabel{ROIidx};
                else
                    tmp = NParticles(:,2:size(Temp.Results.Process.ROIlabel,1)+1);
                    NParticles = [NParticles(:,1) sum(tmp,2)];
                    ROIChoice{end+1,1} = 'All';
                end
            else
                ROIChoice{1} = Temp.Results.Process.ROIlabel{1};
            end
        end
    end
    
    
    
    
    [ImmTracks, Dummy] = calculateImmobileTracks...
        (Tracks, Parameters.ThreshL,...
        Parameters.minBoundFrames, Parameters.ThreshH, 0);
    
    
    if isempty(ImmTracks)
        TrackLengthHist = [0 0];
    else
        
        % Identify bound molecules for each of the movies;
        
        % Calculate survival histogram;
        TrackLengthHist = ...
            calculateTrackLength(ImmTracks, ...
            Parameters.FrameTime,Parameters.minBoundFrames);
        
        
        
    end
    
    
    
else
    errordlg('Some of the files have not been preprocessed');
    return
end

% find the data containing the lowest number of frames;
TimePoints = [];


nMax = length(NParticles(:,1));
TimePoints = NParticles(:,1);


%     if length(NParticles(:,1)) < nMax;
%
%         TimePoints = NParticles(:,1);
%         nMax = length(TimePoints);
%     end


if isempty(TimePoints)
    errordlg('Too few bound particles to produce an histogram');
    return
end

% Accumulate the histogram
Hist_Matrix = zeros(nMax,1);

% for i = 1:nFiles
Hist_Matrix(1:nMax,1) = NParticles(1:nMax,2);
% end

CumNParticles(:,1) = TimePoints * Parameters.FrameTime;
CumNParticles(:,2) = sum(Hist_Matrix, 2);

%Calculate the avg intensity


% Accumulate the different survival histograms;

% find the histogram containing the longest track
nMax = 1;
TimePoints = [];
% for i = 1:nFiles
%
%     if length(TrackLengthHist{i}(:,1)) > nMax;

TimePoints = TrackLengthHist(:,1);
%         nMax = length(TimePoints);
%     end
% end

if isempty(TimePoints)
    errordlg('Too few bound particles to produce an histogram');
    return
end


% Accumulate the histogram
Hist_Matrix = zeros(nMax,1);

% for i = 1:nFiles
lengthHist = length(TrackLengthHist(:,1));
Hist_Matrix(1:lengthHist,1) = TrackLengthHist(:,2);
% end

CumTrackLengthHist(:,1) = TimePoints;
CumTrackLengthHist(:,2) = sum(Hist_Matrix, 2);

% Copy Histogram to handles
handles.CumHist = CumTrackLengthHist;


% Calculate the fraction of bound molecules.

% First calculate the residence time histogram
ResTimeHist(:,1)= CumTrackLengthHist(:,1);
ResTimeHist(:,2) = ...
    [(CumTrackLengthHist(1:end-1,2) - CumTrackLengthHist(2:end,2));CumTrackLengthHist(end,2)];

% Then calculate the total number of bound spots
TotalBoundMolecules = sum(ResTimeHist(:,2).*ResTimeHist(:,1)/Parameters.FrameTime);

% and the total number of spots
TotalMolecules = sum(CumNParticles(:,2));

% Finally divide the two.
PartialBoundFraction = TotalBoundMolecules/TotalMolecules;

BFerror = (sqrt(TotalBoundMolecules)/TotalBoundMolecules + ...
    sqrt(TotalMolecules)/TotalMolecules)*PartialBoundFraction;


% NOTE: THIS IS ONLY THE PARTIAL BOUND FRACTION
% BECAUSE WE ARE LOOKING ONLY AT PARTICLES BOUND
% FOR MORE THAN Nmin. We HAVE TO FIT THE HISTOGRAM OF RESIDENCE
% TIMES AND EXTRAPOLATE THE TRUE VALUE OF THE BOUND FRACTION

% Photobleaching correction


if BleachCurveFromInt == 1
    [PB_curve,hbound,vbound] = PBsections(Temp.Results.Data.imageStack,Temp.Results.Process.ROIimage,1,1,1);
    PB_curve = reshape(PB_curve,numel(PB_curve),1,1);
    CumNParticles(:,2) = PB_curve;
end
if ~isnan(CumNParticles(1,2))
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
    % % PhotobAxes = findobj ('Tag', 'photob_axes');
    % % axes(PhotobAxes);
    %
    %
    % xlabel('Time [s]');
    % ylabel('Normalized Counts');
    %
    %
    %
    % plot(CumNParticles(:,1), CumNParticles(:,2)/max(CumNParticles(:,4)), '+', 'MarkerEdgeColor', [0.5 0.5 0.5]');
    % hold on
    % plot(CumTrackLengthHist(:,1), CumTrackLengthHist(:,2)/max(CumTrackLengthHist(:,2)),'ok');
    % plot(CumNParticles(:,3), CumNParticles(:,4)/max(CumNParticles(:,4)),'r');
    % hold off;
    %
    % box on;
    %
    % title({'Photobleaching kinetics:', ['k_{b1} = ', num2str(BleachRates(1), 3),...
    %     ' s^{-1}, k_{b2} = ', num2str(BleachRates(2), 3), ' s^{-1}, f_1 = ',  ...
    %     num2str(BleachRates(3), 3)]})
    %
    % hold off;
    % legend('Photobleaching Decay', 'Bound molecules Decay');
    
    
    % Copy CumNparticles to handles
    handles.NParticles = CumNParticles;
    handles.BleachRates = BleachRates;
    
    
    % Ask if you want to correct the histogram for photobleaching;
    % answer = questdlg('Do you want to correct the histogram for photobleaching?');
    answer = 'Yes';
    switch answer
        case 'Yes'
            % Recalculate the survival probability corrected for photobleaching
            CumTrackLengthHist(:,2) = CumTrackLengthHist(:,2)./ ...
                (BleachRates(3)*exp(-BleachRates(1).* CumTrackLengthHist(:,1)) + ...
                (1-BleachRates(3))*exp(-BleachRates(2).* CumTrackLengthHist(:,1)));
            % Recalculate the residence time distribution corrected for
            % photobleaching
            ResTimeHist(:,1)= CumTrackLengthHist(:,1);
            ResTimeHist(:,2) = ...
                [(CumTrackLengthHist(1:end-1,2) - CumTrackLengthHist(2:end,2));CumTrackLengthHist(end,2)];
            
            
    end
    % 
    
    
    
    % If the user wanted to analyze the residence time, convert the cumulative
    % histogram to a residence time histogram.
    
    if Parameters.Flag == 1;
        if size(ResTimeHist,1) > 4
            deltaT = ResTimeHist(2,1) - ResTimeHist(1,1);
            % Normalize the Survival probability for the partial bound fraction
            CumTrackLengthHist(:,2) = CumTrackLengthHist(:,2)/CumTrackLengthHist(1,2)...
                .*PartialBoundFraction;
            
            % Normalize the residence time histogram for the partial bound fraction
            ResTimeHist(:,2) = ResTimeHist(:,2)/(deltaT*sum(ResTimeHist(:,2)))*...
                PartialBoundFraction;
            
            
            
            
            
            % Bin the histogram according to the selected settings;
            
            binFactor = Parameters.bin;
            if binFactor == 0
                binFactor = 1;
            end
            nBins = floor(length(ResTimeHist(:,1))/binFactor);
            
            TimeBins = reshape(ResTimeHist(1:nBins*binFactor,1), binFactor, nBins);
            CountsBins = reshape(ResTimeHist(1:nBins*binFactor,2), binFactor, nBins);
            
            if binFactor ~= 1
                resTimeHist_Binned = [];
                resTimeHist_Binned(:,1) = mean(TimeBins);
                resTimeHist_Binned(:,2) = sum(CountsBins);
            else
                resTimeHist_Binned = ResTimeHist;
            end
            
            % Remove values < 0 from the histogram
            idx = find(resTimeHist_Binned(:,2) <= 0);
            resTimeHist_Binned(idx,:) = [];
            if size(resTimeHist_Binned,1) < 4
                Hist = [];
                FitPar = [];
                return
            end
            % Fit the binned residence time histogram with exponentials;
            [fitpar,espSigma, fit]= ExpDecay_fit_resTime(resTimeHist_Binned, 1 );
            [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit_resTime(resTimeHist_Binned, [1 0.1]);
            
            
            
            %
            %
            %     HistAxes = findobj ('Tag', 'hist_axes');
            %     axes(HistAxes);
            %     bar(resTimeHist_Binned(:,1), resTimeHist_Binned(:,2),'EdgeColor', [0 0 0], ...
            %         'FaceColor','none', 'BarWidth', 1, 'LineWidth', 0.75);
            %     hold on;
            %     plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
            %     plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
            %     plot(fit2(:,1),fit2(:,3),'b','LineWidth',1);
            %     plot(fit2(:,1),fit2(:,4),'b--','LineWidth',1);
            %     xlabel('Residence Time [s]')
            %     ylabel('Counts');
            %     title({['Fit of the residence time histogram',...
            %         ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction,3), ' \pm ',...
            %         num2str(BFerror,2)],...
            %         ['One Component fit k_{off} = ',num2str(fitpar(1),3),...
            %         ' \pm ', num2str(espSigma(1),2), ' s; C_{eq} = ', num2str(fitpar(2),2)],...
            %         ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
            %         ' \pm ', num2str(espSigma2(1),2), 's'], ...
            %         ['k_{off, 2} = ', num2str(fitpar2(2),3),...
            %         ' \pm ', num2str(espSigma2(2),2), 's; f_1 = ', num2str(fitpar2(3),3),...
            %         '\pm', num2str(espSigma2(3),2), ' C_{eq} = ',  num2str(fitpar2(4),2)]});
            %     legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit','Component 1', 'Component 2');
            %     box on;
            
            handles.Hist(:, 1:2) = resTimeHist_Binned;
            handles.Hist(:,3) = fit(:,2);
            handles.Hist(:,4) = fit2(:,2);
            handles.Hist(:,5) = fit2(:,3);
            handles.Hist(:,6) = fit2(:,4);
            
            handles.FitPar = [fitpar, fitpar2, PartialBoundFraction; espSigma, espSigma2, BFerror];
        else
            handles.Hist = [];
            handles.FitPar = [];
        end
        
        
        
        
    end
    
    if Parameters.Flag == 2;
        % Logarithmic sampling of the survival prob
        if Parameters.bin ~= 0
            
            
            TimePoints = logspace(log10(min(CumTrackLengthHist(:,1))),...
                log10(max(CumTrackLengthHist(:,1))),Parameters.bin);
            TrLog =[];
            
            TimePoints = round(TimePoints*1000)/1000;
            for i = 1:length(TimePoints)
                
                idx = find(CumTrackLengthHist(:,1) <= TimePoints(i),1,'last');
                TrLog(i,:) = CumTrackLengthHist(idx,:);
            end
            
        else
            TrLog = CumTrackLengthHist;
        end
        
        tpoint1 = TrLog(1,1);
        TrLog(:,1) = TrLog(:,1) - tpoint1;
        
        if size(TrLog,1) > 4
            
            
            % Fit exponentials to the Survival probability
            [fitpar,espSigma, fit]= ExpDecay_fit(TrLog, 1/(tpoint1*10));
            [fitpar2,espSigma2, fit2]= ExpDecay_2Cmp_fit(TrLog, [1/(tpoint1*10) 1/(tpoint1*100)]);
            fit(:,1) = fit(:,1) + tpoint1;
            fit2(:,1) = fit2(:,1) + tpoint1;
            TrLog(:,1) = TrLog(:,1) + tpoint1;
            %     HistAxes = findobj ('Tag', 'hist_axes');
            %     axes(HistAxes);
            %     plot(TrLog(:,1), TrLog(:,2),'ok');
            %     hold on;
            %     plot(fit(:,1),fit(:,2),'g', 'LineWidth', 1);
            %     plot(fit2(:,1),fit2(:,2),'r', 'LineWidth', 1);
            %     plot(fit2(:,1),fit2(:,3),'b', 'LineWidth', 1);
            %     plot(fit2(:,1),fit2(:,4),'b--', 'LineWidth', 1);
            %     xlabel('Time [s]');
            %     ylabel('Counts');
            %     title({['Fit of the survival probability',...
            %         ' - Measured [partial] C_{eq} = ' num2str(PartialBoundFraction,3), ' \pm ',...
            %         num2str(BFerror,2)],...
            %         ['One Component fit k_{off} = ',num2str(fitpar(1),3),...
            %         ' \pm ', num2str(espSigma(1),2), ' s; C_{eq} = ', num2str(fitpar(2),2)],...
            %         ['Two Components fit k_{off,1} = ', num2str(fitpar2(1),3),...
            %         ' \pm ', num2str(espSigma2(1),2), 's'], ...
            %         ['k_{off, 2} = ', num2str(fitpar2(2),3),...
            %         ' \pm ', num2str(espSigma2(2),2), 's; f_1 = ', num2str(fitpar2(3),3),...
            %         '\pm', num2str(espSigma2(3),2), 'C_{eq} = ',  num2str(fitpar2(4),2)]});
            %     legend('Residence Times', 'Single exponential fit', 'Double Exponential Fit');
            %     box on;
            handles.Hist(:, 1:2) = TrLog;
            handles.Hist(:,3) = fit(:,2);
            handles.Hist(:,4) = fit2(:,2);
            handles.Hist(:,5) = fit2(:,3);
            handles.Hist(:,6) = fit2(:,4);
            handles.FitPar = [fitpar, fitpar2, PartialBoundFraction; espSigma, espSigma2, BFerror];
        else
            handles.Hist = [];
            handles.FitPar = [];
        end
    end
    
else
    handles.Hist = [];
    handles.FitPar = [];
end
Hist = handles.Hist;
FitPar = handles.FitPar;

