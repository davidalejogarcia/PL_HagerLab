function [peak_summary_new peak_fit x_sub y_sub ] = peak_seqFit_psf(stack,peak_summary,dx,dy,flagGauss,flagOutput)

% flagGauss = 1 ... pixelated Gaussian
% flagGauss = 2 ... continous Gaussian

if nargin <=2
    dx = 5;
    dy = 5;
end


[nPeaks, nCol] = size(peak_summary);
peak_summary_new = zeros(nPeaks,nCol+6);

peak_summary_new(:,1:nCol) = peak_summary;

peak_fit = struct('raw', {}, 'fit', {});

for iPeak = 1:nPeaks
    
    disp(['Peak dectection - fit with PSF: peak # ',num2str(iPeak),' of ',num2str(nPeaks)])
      
    iStack = peak_summary(iPeak,6);
    img    = double(stack(iStack).data);

    x_center = peak_summary(iPeak,1);
    y_center = peak_summary(iPeak,2);
    
    
     if x_center ~= 0 && y_center ~= 0
        
        %- Extract subimage
        x_center_round = round(x_center);
        y_center_round = round(y_center);
  
        y_sub = (y_center_round-dx:y_center_round+dx);
        x_sub = (x_center_round-dy:x_center_round+dy);
        img_sub = img(y_sub,x_sub);
        
        % Refine x-coordinate of the peak
        x_proj = mean(img_sub,1);  % calculate the x-projection of the sub-image
        Coef_0 = [max(x_proj), 2, x_center, min(x_proj)]; %set initial values for the fit
        Coef_x = nlinfit(x_sub,x_proj,@gaussian_1D_fun, Coef_0); % fit x - projection
        x_out = gaussian_1D_fun(Coef_x,x_sub);
        
        % Refine y-coordinate of the peak
        y_proj = mean(img_sub,2)';  % calculate the y-projection of the sub-image
        Coef_0 = [max(y_proj), 2, y_center, min(y_proj)]; %set initial values for the fit
        Coef_y = nlinfit(y_sub,y_proj,@gaussian_1D_fun, Coef_0); % fit y - projection
        y_out = gaussian_1D_fun(Coef_y,y_sub);
        
        % Fit the 2-D sub image for maxIntensity, sigma and bkground 
        x_center_fit = Coef_x(3);
        y_center_fit = Coef_y(3);
        sigmaStart = (Coef_x(2)+Coef_y(2))/2;
        bgdStart =   (Coef_x(4)+Coef_y(4))/2;
        
        parMod_1    =  {1, x_sub, y_sub, x_center_fit, y_center_fit};
        parStart_1 =  [max(max(img)),sigmaStart,bgdStart];
    
        if      flagGauss == 1
            [parFit_1 ssr_1 psf_fit_1] = gauss_step_fit(img_sub,parMod_1,parStart_1);
        elseif  flagGauss == 2
            [parFit_1 ssr_1 psf_fit_1] = gauss_2D_cont_fit(img_sub,parMod_1,parStart_1);
        end
        
        peak_summary_new(iPeak,nCol+1:end) = [parFit_1 x_center_fit y_center_fit ssr_1];    
        peak_fit(iPeak).raw = img_sub;
        peak_fit(iPeak).fit = psf_fit_1;
        
        
      if flagOutput
        
        figure;
        subplot(2,1,1);
        
        plot (x_sub, x_proj, 'ok');
        hold on;
        plot (x_sub,x_out,'r');
        
        subplot(2,1,2)
        plot (y_sub, y_proj, 'ok');
        hold on;
        plot (y_sub,y_out,'r');
        end
        
        
    end
        
        
        
        
    
end