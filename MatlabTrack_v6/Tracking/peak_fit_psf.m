function [peak_summary_new] = peak_fit_psf(stack,peak_summary,dx,dy)



if nargin <=2
    dx = 5;
    dy = 5;
end



% Find Number of Images
nFrames =  max(peak_summary(:,6));
% initialize the meshgrids for the calculation of gaussians.
sizeFrame = size(stack(1).data);
[x,y] = meshgrid(1:sizeFrame(2), 1:sizeFrame(1));
peak_summary_new =[];

% Size of image

% Cycle over the images;

for iFrame = 1:nFrames;
    % Find peaks in current frame
    disp(['Detecting Peaks in frame ', num2str(iFrame)]);
    idxPeaks = peak_summary(:, 6) == iFrame;   
    nPeaks = sum(idxPeaks);
    if nPeaks > 0
        PeaksInCurrFrame = sortrows(peak_summary(idxPeaks, :), -5);
        PeaksInCurrFrame = [PeaksInCurrFrame(:,1:7), zeros(nPeaks,6)];
        
        
        % Detect if any particle is there at the nth Frames
        if PeaksInCurrFrame(1,1) == 0 && PeaksInCurrFrame(1,2) == 0;
            disp(['No Peaks in frame ', num2str(iFrame)]);
            disp('________________________');
            
        else
            % get image
            img  = double(stack(iFrame).data);
            
            disp([num2str(nPeaks),' Peaks detected in frame ', num2str(iFrame)]);
            
            for iPeak = 1:nPeaks
                % Fit gaussian to the peak
                x_center = PeaksInCurrFrame(iPeak,1);
                y_center = PeaksInCurrFrame(iPeak,2);
                
                
                if x_center ~= 0
                    
                    
                    x_center_round = round(x_center);
                    y_center_round = round(y_center);
                    y_sub = (y_center_round-dy:y_center_round+dy);
                    x_sub = (x_center_round-dx:x_center_round+dx);
                    
                    img_sub = img(y_sub,x_sub);
                    coordinates = {x_sub, y_sub};
                    disp(['Fitting Gaussian to peak ', num2str(iPeak), ' in frame ', num2str(iFrame)]);
                    
                    [parFit, ssr] = gauss_fit2D_cov(img_sub,coordinates);
                    
                    
                    % Add fit parameters to the output
                    PeaksInCurrFrame(iPeak, 8:13) = [parFit, ssr];
                    %                 % Test plot
                    %                 figure;
                    %                 subplot(2,1,1);
                    %                 imagesc(img);
                    %                 hold on;
                    %                 colormap(gray);
                    %                 plot(PeaksInCurrFrame(10, iPeak),PeaksInCurrFrame(11, iPeak), 'or');
                    %
                    % Subtract found particle from image;
                    exponent = ((x - parFit(4))/parFit(2)).^2 +...
                        ((y - parFit(5))/parFit(2)).^2;
                    
                    img_del = parFit(1) * exp(-exponent/2);
                    
                    img =img - img_del;
                    
                    %                 % Test plot
                    %                 subplot(2,1,2);
                    %                 imagesc(img);
                    %                 hold on;
                    %                 colormap(gray);
                    %                 plot(PeaksInCurrFrame(10, iPeak),PeaksInCurrFrame(11, iPeak), 'or');
                    %                 hold off;
                    
                end
                
            end
            
        end
        
        peak_summary_new = [peak_summary_new; PeaksInCurrFrame];
    else
        disp(['No Peaks in frame ', num2str(iFrame)]);
        disp('________________________');
    end
   
   
     
end
%move ROI data to last column to preserve downstream processing

peak_summary_new = [peak_summary_new(:,1:6) peak_summary_new(:,8:13) peak_summary_new(:,7)];
peak_summary_new = sortrows(peak_summary_new,[6 13]);
TotPart = sum(peak_summary_new(:,1) ~= 0);
disp([num2str(TotPart), ' Particles found in ', num2str(nFrames), ' Frames']);

%% OLD CODE
% 
% [nPeaks, nCol] = size(peak_summary);
% peak_summary_new = zeros(nPeaks,nCol+6);
% 
% peak_summary_new(:,1:nCol) = peak_summary;   
%     

% for iPeak = 1:nPeaks
%     
%     disp(['Peak dectection - fit with PSF: peak # ',num2str(iPeak),' of ',num2str(nPeaks)])
%       
%     iStack = peak_summary(iPeak,6);
%     img    = double(stack(iStack).data);
% 
%     x_center = peak_summary(iPeak,1);
%     y_center = peak_summary(iPeak,2);
%  
%     if x_center ~= 0 && y_center ~= 0
%         
%     %- Extract subimage
%     x_center_round = round(x_center);
%     y_center_round = round(y_center);
%   
%     y_sub = (y_center_round-dy:y_center_round+dy);
%     x_sub = (x_center_round-dx:x_center_round+dx);
%     img_sub = img(y_sub,x_sub);
%     
%     
%     % My fitting
%     coordinates = {x_sub, y_sub};
%     
%     [parFit, ssr] = gauss_fit2D_cov(img_sub,coordinates);
%     
%     
% 
%     peak_summary_new(iPeak,nCol+1:end) = [parFit ssr];    
% 
%     end
% end

