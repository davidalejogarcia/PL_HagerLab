function [peak_summary_new peak_fit x_sub y_sub ] = peak_Fit_psf_noXY(stack,peak_summary,dx,dy,flagGauss,flagOutput)

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
        
     

        
        % Fit the 2-D sub image for maxIntensity, sigma and bkground 

        sigmaStart =  2;
        bgdStart =   (mean(img_sub(:,1))+mean(img_sub(:,end))+mean(img_sub(1,:))+mean(img_sub(end,:)))/4;
        
        parMod_1    =  {1, x_sub, y_sub, x_center, y_center};
        parStart_1 =  [max(max(img)),sigmaStart,bgdStart];
    
        if      flagGauss == 1
            [parFit_1 ssr_1 psf_fit_1] = gauss_step_fit(img_sub,parMod_1,parStart_1);
        elseif  flagGauss == 2
            [parFit_1 ssr_1 psf_fit_1] = gauss_2D_cont_fit(img_sub,parMod_1,parStart_1);
        end
        
        peak_summary_new(iPeak,nCol+1:end) = [parFit_1 x_center y_center ssr_1];    
        peak_fit(iPeak).raw = img_sub;
        peak_fit(iPeak).fit = psf_fit_1;
        
        
        if flagOutput

        %- Plot results
            figure
            subplot(2,1,1)
            hold on
            h=surf(double(img_sub));
            v= axis;
            set(h,'EdgeColor','k');
            set(gca,'XGrid','off')
            set(gca,'YGrid','off')
            set(gca,'ZGrid','off')
     
            subplot(2,1,2)     
            h=surf(double(psf_fit_1));
            set(h,'EdgeColor','k');
            set(gca,'XGrid','off')
            set(gca,'YGrid','off')
            set(gca,'ZGrid','off')
            axis(v)
        end



        
    end
        
        
        
        
    
end