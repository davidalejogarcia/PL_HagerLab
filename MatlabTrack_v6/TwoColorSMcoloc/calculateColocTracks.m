function [boundParticles, allboundParticles] = calculateColocTracks(T1,T2,ThreshL,minBound)

D_mat = [];
allboundParticles = [];
boundParticles = [];
%obtain all points that are within ThreshL
SegmentID = 0;
nImages = max(max(T1(:,3)),max(T2(:,3)));

for i = 1:nImages
    T1_tmp = T1(T1(:,3) == i,:);
    T2_tmp = T2(T2(:,3) == i,:);
    for j = 1:size(T1_tmp,1)
        T1_1 = T1_tmp(j,:);
        for k = 1:size(T2_tmp,1)
            T2_1 = T2_tmp(k,:);
            D_cur = sqrt(((T1_1(1) - T2_1(1)).^2) + ((T1_1(2) - T2_1(2)).^2));
            if D_cur < ThreshL
                D_cur1 = [i, T1_1(4), T2_1(4), D_cur];
                D_mat = [D_mat; D_cur1];
            end
        end
    end
end
if ~isempty(D_mat)

    p1 = unique(D_mat(:,2));
    
    for i = 1:length(p1)
        cur1 = D_mat(D_mat(:,2) == p1(i),:);
        p2 = unique(cur1(:,3));
        for j = 1:length(p2)
            cur2 = cur1(cur1(:,3) == p2(j),:);
            if size(cur2,1) > 1
                Tstep = diff(cur2(:,1));
                jumps = find(Tstep > 1);
                if ~isempty(jumps)
                    for k = 1:length(jumps)
                        if k == 1
                            ind1 = 1;
                        else
                            ind1 = jumps(k)+1;
                        end
                        ind2 = jumps(k);
                        SegmentID = SegmentID + 1;

                        catTrack = [cur2(ind1:ind2,1:4),SegmentID*ones(length(ind1:ind2),1)];
                        allboundParticles = [allboundParticles; catTrack];
                    end
                else
                    SegmentID = SegmentID + 1;
                    catTrack = [cur2(:,1:4),SegmentID*ones(size(cur2,1),1)];
                    allboundParticles = [allboundParticles; catTrack];
                end
            else
                SegmentID = SegmentID + 1;
                catTrack = [cur2(:,1:4),SegmentID];
                allboundParticles = [allboundParticles; catTrack];
            end
        end
    end

    %Discard segments that are shorter than minBound
    SavedSegID = 0;
    for i = 1:max(allboundParticles(:,5))
        tmp = allboundParticles(allboundParticles(:,5) == i,:);
        if size(tmp,1) >= minBound
            SavedSegID = SavedSegID + 1;
            curTrack = [tmp(:,1:4), SavedSegID*ones(size(tmp,1),1)];
            boundParticles = [boundParticles; curTrack];
        end
    end
end
                
                
            