function BatchTrack(lpass,hpass,threshold,windowSz,isFitPSF,maxJump,shTrack,closeGaps)

[fnames,pnames] = uigetfile('*.tif','Select the images to open','','MultiSelect','on');
savefold = uigetdir(pwd,'Select a save location');

if ~iscell(fnames)
    fnames = {fnames};
end
if ~iscell(pnames)
    pnames = {pnames};
end
if ~isempty(fnames)
    for i = 1:length(fnames)
        clear Results;
        Results.Data.fileName = fnames{i};
        Results.Data.pathName = pnames{1};
        
        %Load in the image
        [Results.Data.imageStack,Results.Data.nImages] = TIFread([pnames{1}, fnames{i}]);
        
        Results.Data.clims = [min(min(Results.Data.imageStack(1).data))...
            max(max(Results.Data.imageStack(1).data))];
        
        %filter the image
        
        for j =  1:Results.Data.nImages
            Results.Process.filterStack(j).data = ...
                bpass(Results.Data.imageStack(j).data, lpass, hpass);
        end
        Results.Process.clims = [min(min(Results.Process.filterStack(1).data))...
            max(max(Results.Process.filterStack(1).data))];
        %Generate the ROI
        
        rmax = 1;
        
        pxSize = 0.104;
        
        Results.Process.ROIlabel{1,:} = 'ROI1';
        centX = Results.Data.imageStack(1).width/2;
        centY = Results.Data.imageStack(1).height/2;
        rmax_px = (rmax/pxSize) + 1 + 3;
        xdata = (centX - rmax_px):0.1:(centX +rmax_px);
        ydata1 = centY + sqrt((rmax_px^2 - (xdata - centX).^2));
        ydata2 = centY - sqrt((rmax_px^2 - (xdata - centX).^2));
        Results.Process.ROIpos{1,:} = [xdata',ydata1';flipud(xdata'), ydata2';xdata(1),ydata1(1)];
        
        Results.Process.ROIimage{1,:} = false(size(Results.Data.imageStack(1).data));
        for m = floor(centX - rmax_px):ceil(centX + rmax_px)
            for n = floor(centY - rmax_px):ceil(centY + rmax_px)
                if sqrt((centX - m).^2 + (centY - n).^2) <= rmax_px
                    Results.Process.ROIimage{1,:}(m,n) = 1;
                end
            end
        end
        Results.Parameters.Tracking = [lpass,hpass,threshold,windowSz,maxJump,closeGaps,shTrack];
        Results.Parameters.Acquisition.pixelSize = 0.115;
        Results.Parameters.Acquisition.frameTime = 0.2;
        
        %Find Particles
        Results.Tracking.Centroids = findParticles(Results.Process.filterStack, threshold, hpass,windowSz);
        if max(Results.Tracking.Centroids(:,1)) > 0
            if isFitPSF
                CentroidInRoi = InsideROIcheck2(Results.Tracking.Centroids, Results.Process.ROIimage);

                Centroid = CentroidInRoi;
                Particles = peak_fit_psf(Results.Data.imageStack,...
                    Centroid,windowSz,windowSz);
                Particles2 = InsideROIcheck2(Particles,Results.Process.ROIimage);
                Results.Tracking.Particles = Particles2;
            end
            Results.Tracking.Peaks = peaks;

            %Track 
            Trackparam.mem = closeGaps;
            Trackparam.good = shTrack;
            Trackparam.dim         =  2;
            Trackparam.quiet       =  0;
            if isFitPSF % if particle position has been evaluated via PSF fitting
                Particles = Results.Tracking.Particles(:,[10 11 6 13]);
            else
                Particles = Results.Tracking.Centroids(:,[1 2 6 7]);
            end;
            [Tracks, TrkPtsAdded, errorcode] = trackfunctIG(Particles(:,1:3),maxJump,Trackparam);
            if errorcode == 0 && ~isempty(Tracks)
                Tracks = InsideROIcheck2(Tracks,Results.Process.ROIimage);
                Results.Tracking.Tracks = Tracks;
                if isFitPSF
                    ParticlesNew = Results.Tracking.Particles;

                    x_ind = 10;
                    y_ind = 11;
                else
                    ParticlesNew = Results.Tracking.Centroids;
                    x_ind = 1;
                    y_ind = 2;
                end
                Particles = ParticlesNew;

                for j = 1:size(Tracks,1)
                    x_pos = Tracks(j,1);
                    y_pos = Tracks(j,2);
                    frame_num = Tracks(j,3);

                    pIx1 = find(Particles(:,x_ind) == x_pos & ...
                        Particles(:,y_ind) == y_pos & ...
                        Particles (:,6) == frame_num);
                    if isempty(pIx1)
                        ParticleAdd(:,1) = x_pos;
                        ParticleAdd(:,2) = y_pos;
                        ParticleAdd(:,6) = frame_num;
                        if isFitPSF
                            ParticleAdd(:,10:11) = ParticleAdd(:,1:2);
                            ParticleAdd(:,12) = 0;
                            if isfield(Results.Process,'ROIpos')
                                ParticleAdd(:,13) = 0;
                            end
                        else
                            if isfield(Results.Process,'ROIpos')
                                ParticleAdd(:,7) = 0;
                            end
                        end
                        ParticlesNew = [ParticlesNew; ParticleAdd];
                    end


                end
                ParticlesNew = InsideROIcheck2(ParticlesNew,Results.Process.ROIimage);
                if isFitPSF
                    ParticlesNew = sortrows(ParticlesNew,[6,13]);
                    Results.Tracking.Particles = ParticlesNew;
                else
                    ParticlesNew = sortrows(ParticlesNew,[6,7]);
                    Results.Tracking.Centroids = ParticlesNew;
                end


                %PreAnalysis
                Tracks = Results.Tracking.Tracks;
                if isFitPSF
                    Particles = Results.Tracking.Particles;
                else
                    Particles = Results.Tracking.Centroids;
                end

                [Results.PreAnalysis.Tracks_um, Results.PreAnalysis.NParticles, Results.PreAnalysis.IntensityHist] = preProcess_noGUI(Tracks,Results.Data.imageStack,...
                    Particles, pxSize, Results.Data.nImages, Results.Data.fileName,Results.Process.ROIpos);
                Results.isFitPSF = isFitPSF;
                Results.Analysis = [];
                Version = 2;
                save([savefold, filesep,fnames{i}(1:end-4),'_tracked_preprocess.mat'],'Results','Version');
            end
        end
        
    end
end