function [closeParticles, allboundParticles] = calculateColocParticles(P1,P2,ThreshL)

closeParticles = [];
allboundParticles = [];
%obtain all points that are within ThreshL
% SegmentID = 0;
nImages = max(max(P1(:,6)),max(P2(:,6)));
if size(P1,2) > 8 && size(P2,2) > 8
    x_ind = 10;
    y_ind = 11;
else
    x_ind = 1;
    y_ind = 2;
end
for i = 1:nImages
    P1_tmp = P1(P1(:,6) == i,:);
    P2_tmp = P2(P2(:,6) == i,:);
    for j = 1:size(P1_tmp,1)
        P1_1 = P1_tmp(j,:);
        for k = 1:size(P2_tmp,1)
            P2_1 = P2_tmp(k,:);
            D_cur = sqrt(((P1_1(x_ind) - P2_1(x_ind)).^2) + ((P1_1(y_ind) - P2_1(y_ind)).^2));
            if D_cur < ThreshL
                D_cur1 = [i, j, k, D_cur];
                closeParticles = [closeParticles; D_cur1];
            end
        end
    end
end
% 
% p1 = unique(D_mat(:,2));
% boundParticles = [];
% for i = 1:length(p1)
%     cur1 = D_mat(D_mat(:,2) == p1(i),:);
%     p2 = unique(cur1(:,3));
%     for j = 1:length(p2)
%         cur2 = cur1(cur1(:,3) == p2(j),:);
%         if size(cur2,1) > 1
%             Tstep = diff(cur2(:,1));
%             jumps = find(Tstep > 1);
%             if ~isempty(jumps)
%                 for k = 1:length(jumps)
%                     if k == 1
%                         ind1 = 1;
%                     else
%                         ind1 = jumps(k)+1;
%                     end
%                     ind2 = jumps(k);
%                     SegmentID = SegmentID + 1;
% 
%                     catTrack = [cur2(ind1:ind2,1:4),SegmentID*ones(length(ind1:ind2),1)];
%                     allboundParticles = [allboundParticles; catTrack];
%                 end
%             else
%                 SegmentID = SegmentID + 1;
%                 catTrack = [cur2(:,1:4),SegmentID*ones(size(cur2,1),1)];
%                 allboundParticles = [allboundParticles; catTrack];
%             end
%         else
%             SegmentID = SegmentID + 1;
%             catTrack = [cur2(:,1:4),SegmentID];
%             allboundParticles = [allboundParticles; catTrack];
%         end
%     end
% end
% 
% %Discard segments that are shorter than minBound
% SavedSegID = 0;
% for i = 1:max(allboundParticles(:,5))
%     tmp = allboundParticles(allboundParticles(:,5) == i,:);
%     if size(tmp,1) >= minBound
%         SavedSegID = SavedSegID + 1;
%         curTrack = [tmp(:,1:4), SavedSegID*ones(size(tmp,1),1)];
%         boundParticles = [boundParticles; curTrack];
%     end
% end
                
                
            