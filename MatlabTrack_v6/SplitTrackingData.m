function Results_out = SplitTrackingData(Results,hor_bound,ver_bound)

Results_out = cell(size(ver_bound,1),size(hor_bound,1));

for i = 1:size(ver_bound,1)
    for j = 1:size(hor_bound,1)
        Results_out{i,j} = Results;
        for k = 1:Results.Data.nImages
            Results_out{i,j}.Data.imageStack(k).data = Results.Data.imageStack(k).data(ver_bound(i,1):ver_bound(i,2),hor_bound(j,1):hor_bound(j,2));
            Results_out{i,j}.Data.imageStack(k).height = size(Results_out{i,j}.Data.imageStack(k).data,1);
            Results_out{i,j}.Data.imageStack(k).width = size(Results_out{i,j}.Data.imageStack(k).data,2);
            Results_out{i,j}.Process.filterStack(k).data = Results.Process.filterStack(k).data(ver_bound(i,1):ver_bound(i,2),hor_bound(j,1):hor_bound(j,2));
            
        end
        if isfield(Results.Process,'ROIimage')
            for k = 1:length(Results.Process.ROIimage)
                Results_out{i,j}.Process.ROIimage{k,:} = Results.Process.ROIimage{k}(ver_bound(i,1):ver_bound(i,2),hor_bound(j,1):hor_bound(j,2));
                Results_out{i,j}.Process.ROIpos{k,:}(:,1) = Results.Process.ROIpos{k}(:,1) - hor_bound(j,1) + 1;
                Results_out{i,j}.Process.ROIpos{k,:}(:,2) = Results.Process.ROIpos{k}(:,2) - ver_bound(i,1) + 1;
                Results_out{i,j}.Process.ROIpos{k,:} = [min(max(1,Results_out{i,j}.Process.ROIpos{k,:}(:,1)),Results_out{i,j}.Data.imageStack(1).width),...
                    min(max(1,Results_out{i,j}.Process.ROIpos{k,:}(:,2)),Results_out{i,j}.Data.imageStack(1).height)];
            end
        end
        if isfield(Results,'Tracking')
            if isfield(Results.Tracking,'Particles')
                Particles_all = Results.Tracking.Particles;
                ind_keep = find(Particles_all(:,1) >= hor_bound(j,1) & Particles_all(:,1) <= hor_bound(j,2) & Particles_all(:,2) >= ver_bound(i,1) & Particles_all(:,2) <= ver_bound(i,2));
                Results_out{i,j}.Tracking.Particles = Particles_all(ind_keep,:);
                Results_out{i,j}.Tracking.Particles(:,1) = Results_out{i,j}.Tracking.Particles(:,1) - hor_bound(j,1) + 1;
                Results_out{i,j}.Tracking.Particles(:,2) = Results_out{i,j}.Tracking.Particles(:,2) - ver_bound(i,1) + 1;
                Results_out{i,j}.Tracking.Particles(:,10) = Results_out{i,j}.Tracking.Particles(:,10) - hor_bound(j,1) + 1;
                Results_out{i,j}.Tracking.Particles(:,11) = Results_out{i,j}.Tracking.Particles(:,11) - ver_bound(i,1) + 1;
            end
            if isfield(Results.Tracking,'Centroids')
                
                Particles_all = Results.Tracking.Centroids;
                ind_keep = find(Particles_all(:,1) >= hor_bound(j,1) & Particles_all(:,1) <= hor_bound(j,2) & Particles_all(:,2) >= ver_bound(i,1) & Particles_all(:,2) <= ver_bound(i,2));
                Results_out{i,j}.Tracking.Centroids = Particles_all(ind_keep,:);
                Results_out{i,j}.Tracking.Centroids(:,1) = Results_out{i,j}.Tracking.Centroids(:,1) - hor_bound(j,1) + 1;
                Results_out{i,j}.Tracking.Centroids(:,2) = Results_out{i,j}.Tracking.Centroids(:,2) - ver_bound(i,1) + 1;
            end
            if isfield(Results.Tracking,'Tracks')
                Tracks_all = Results.Tracking.Tracks;
                ind_keep = find(Tracks_all(:,1) >= hor_bound(j,1) & Tracks_all(:,1) <= hor_bound(j,2) & Tracks_all(:,2) >= ver_bound(i,1) & Tracks_all(:,2) <= ver_bound(i,2));
                Results_out{i,j}.Tracking.Tracks = Tracks_all(ind_keep,:);
                Results_out{i,j}.Tracking.Tracks(:,1) = Results_out{i,j}.Tracking.Tracks(:,1) - hor_bound(j,1) + 1;
                Results_out{i,j}.Tracking.Tracks(:,2) = Results_out{i,j}.Tracking.Tracks(:,2) - ver_bound(i,1) + 1;
                if isfield(Results,'PreAnalysis')
                    if isfield(Results.PreAnalysis,'Tracks_um')
                        Results_out{i,j}.PreAnalysis.Tracks_um = Results.PreAnalysis.Tracks_um(ind_keep,:);
                        Results_out{i,j}.PreAnalysis.Tracks_um(:,1) = Results_out{i,j}.PreAnalysis.Tracks_um(:,1) - (hor_bound(j,1) + 1)*Results_out{i,j}.Parameters.Acquisition.pixelSize;
                        Results_out{i,j}.PreAnalysis.Tracks_um(:,2) = Results_out{i,j}.PreAnalysis.Tracks_um(:,2) - (ver_bound(i,1) + 1)*Results_out{i,j}.Parameters.Acquisition.pixelSize;
                    end
                    if ~isfield(Results.Tracking,'CheckParticles')
                        for p = 1:max(Results.PreAnalysis.NParticles(:,1))
                            Results_out{i,j}.PreAnalysis.NParticles(p,2) = size(Results_out{i,j}.Tracking.Particles(Results_out{i,j}.Tracking.Particles(:,6) == p,:),1);
                        end
                    end
                end
                
            end
            if isfield(Results.Tracking,'CheckTracks')
                Tracks_all = Results.Tracking.CheckTracks;
                ind_keep = find(Tracks_all(:,1) >= hor_bound(j,1) & Tracks_all(:,1) <= hor_bound(j,2) & Tracks_all(:,2) >= ver_bound(i,1) & Tracks_all(:,2) <= ver_bound(i,2));
                Results_out{i,j}.Tracking.CheckTracks = Tracks_all(ind_keep,:);
                Results_out{i,j}.Tracking.CheckTracks(:,1) = Results_out{i,j}.Tracking.CheckTracks(:,1) - hor_bound(j,1) + 1;
                Results_out{i,j}.Tracking.CheckTracks(:,2) = Results_out{i,j}.Tracking.CheckTracks(:,2) - ver_bound(i,1) + 1;
            end
            if isfield(Results.Tracking,'CheckParticles')
                Particles_all = Results.Tracking.CheckParticles;
                ind_keep = find(Particles_all(:,1) >= hor_bound(j,1) & Particles_all(:,1) <= hor_bound(j,2) & Particles_all(:,2) >= ver_bound(i,1) & Particles_all(:,2) <= ver_bound(i,2));
                Results_out{i,j}.Tracking.CheckParticles = Particles_all(ind_keep,:);
                Results_out{i,j}.Tracking.CheckParticles(:,1) = Results_out{i,j}.Tracking.CheckParticles(:,1) - hor_bound(j,1) + 1;
                Results_out{i,j}.Tracking.CheckParticles(:,2) = Results_out{i,j}.Tracking.CheckParticles(:,2) - ver_bound(i,1) + 1;
                Results_out{i,j}.Tracking.CheckParticles(:,10) = Results_out{i,j}.Tracking.CheckParticles(:,10) - hor_bound(j,1) + 1;
                Results_out{i,j}.Tracking.CheckParticles(:,11) = Results_out{i,j}.Tracking.CheckParticles(:,11) - ver_bound(i,1) + 1;
                if isfield(Results,'PreAnalysis')
                    for p = 1:max(Results.PreAnalysis.NParticles(:,1))
                        Results_out{i,j}.PreAnalysis.NParticles(p,2) = size(Results_out{i,j}.Tracking.CheckParticles(Results_out{i,j}.Tracking.CheckParticles(:,6) == p,:),1);
                    end
                end
            end
        end
    end
end