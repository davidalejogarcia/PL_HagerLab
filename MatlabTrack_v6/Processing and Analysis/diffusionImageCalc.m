function D = diffusionImageCalc(Results,R2thresh,frameThresh,varargin)

%Calculates the diffusion coefficents of individual tracks by fitting the
%total displacement to a line. The diffusion coefficient is the slope of
%the line. The input R2thresh determines how good the fit for an individual
%track should be in order to keep it. If frameTime is not specified, the
%parameter contained within Results will be used.

if isempty(varargin)
    frameTime = Results.Parameters.Acquisition.frameTime;
else
    frameTime = varargin{1};
end

Tracks = Results.PreAnalysis.Tracks_um;
if length(varargin) == 2
    ROIvec = varargin{2};
    Tracks2 = [];
    for i = 1:length(ROIvec)
        tracks_tmp = Tracks(Tracks(:,5) == ROIvec(i),:);
        Tracks2 = [Tracks2; tracks_tmp];
    end
    Tracks = Tracks2;
end

D = zeros(max(Tracks(:,4)),1);
for i = 1:max(Tracks(:,4))
    
    track1 = Tracks(Tracks(:,4) == i,:);
    if size(track1,1) > frameThresh
        dxy_t = zeros(size(track1,1),2);
        dxy_t(:,1) = track1(:,1) - track1(1,1);
        dxy_t(:,2) = track1(:,2) - track1(1,2);
        dxy2_t = dxy_t.^2;
        disp = sqrt(sum(dxy2_t,2));
        t = track1(:,3);
        t = (t-t(1))*frameTime;
        [L_Coef, L_Sigma, L_Fit] = Line_fit([t,disp], [1 1]);
        dbar = mean(disp);
        SST = sum((disp - dbar).^2);
        SSR = sum((disp - L_Fit(:,2)).^2);
        if SSR/SST <= (1-R2thresh)
            D(i) = L_Coef(1);
        end
    end
end

D(D == 0) = [];
