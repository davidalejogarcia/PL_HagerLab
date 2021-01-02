function Centroid_cull = cullCentroids(CentroidsInRoi, ROI_nest_matrix)

%Takes the input cell array CentroidsInRoi, in which each cell contains the
%data from a single ROI, and generates a single 7-column matrix where each
%Particle is represented once. This is necessary because if a particle is
%inside the inner-most ROI, it will show up in all lists. Uses the
%ROI_nest_matrix which is generated from roiNest.m

Centroids_all = [];
for i = 1:length(CentroidsInRoi)
    Centroids_all = [Centroids_all; CentroidsInRoi{i}];
end

if max(ROI_nest_matrix(:)) > 0
    Centroid_cull = [];
    rem_ind = [];
    for i = 1:length(Centroids_all)-1
            if isempty(rem_ind) || max(i == rem_ind) == 0 
                Centroid_cur = Centroids_all(i,:);
                Centroids_rest = Centroids_all(i+1:end,:);

                Centroid_diff = abs(Centroids_rest - repmat(Centroid_cur,size(Centroids_rest,1),1));

                repeatIndx = find(Centroid_diff(:,1) <= eps & Centroid_diff(:,2) <= eps  ...
                    & Centroid_diff(:,3) <= eps & Centroid_diff(:,4) <= eps ... 
                    & Centroid_diff(:,5) <= eps & Centroid_diff(:,6) <= eps);
                if ~isempty(repeatIndx)
                    Centroid_keep = Centroid_cur(:,1:6);
                    whichROI = [Centroid_cur(7); Centroids_rest(repeatIndx,7)];
                    j = 1;
                    inner = 1;
                    while j <= length(whichROI)
                        test = ROI_nest_matrix(whichROI(j),:);
                        inner_tmp = find(test == 1,1);
                        
                        if ~isempty(inner_tmp) && ~isempty(intersect(inner_tmp,whichROI))
                            j = inner_tmp;
                            inner = j;
                        else
                            j = j+1;
                        end
                    end
                    Centroid_keep(7) = whichROI(inner);
                    Centroid_cull = [Centroid_cull; Centroid_keep];
                    
                else
                    Centroid_cull = [Centroid_cull; Centroid_cur];
                end
                rem_ind = [rem_ind; repeatIndx];
            end
    end
else
    Centroid_cull = Centroids_all;
end

