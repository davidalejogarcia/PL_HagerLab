function [removed, merged] = cullTracksSingle(inTrackDim,outTracks,wholeImage,curInd)



%create a 3D image from the input stack
% for i = 1:length(imStack)
%     wholeImage(:,:,i) = imStack(i).data;
% end

removed = [];
merged = [];
nTracks = length(outTracks);

for i = curInd + 1:nTracks
    checkTrack = outTracks{i};
    
    tvec1 = (inTrackDim(1,5):inTrackDim(1,5)+inTrackDim(1,2))';
    if inTrackDim(1,4) == 1
        xvec1 = (inTrackDim(1,6):inTrackDim(1,6)+inTrackDim(1,3))';
        yvec1 = (inTrackDim(1,7):inTrackDim(1,8))';
    else
        yvec1 = (inTrackDim(1,6):inTrackDim(1,6)+inTrackDim(1,3))';
        xvec1 = (inTrackDim(1,7):inTrackDim(1,8))';
    end
    if (outTracks{i}(1,3) >= tvec1(1) && outTracks{i}(1,3) <= tvec1(end) && ...
            outTracks{i}(1,1) >= xvec1(1) && outTracks{i}(1,1) <= xvec1(end) && ...
            outTracks{i}(1,2) >= yvec1(1) && outTracks{i}(1,2) <= yvec1(end))
        twin = (tvec1(1):tvec1(end))';
        curImage = wholeImage(:,:,twin);
        xwin = (xvec1(1):xvec1(end))';
        ywin = (yvec1(1):yvec1(end))';
        xProj = max(permute(curImage(ywin,xwin,:),[2,3,1]),[],3);
        yProj = max(permute(curImage(ywin,xwin,:),[1,3,2]),[],3);
        
        xbox1 = [tvec1(1),xvec1(1); tvec1(1),xvec1(end); tvec1(end),xvec1(end); tvec1(end),xvec1(1); tvec1(1), xvec1(1)];
        ybox1 = [tvec1(1),yvec1(1); tvec1(1),yvec1(end); tvec1(end),yvec1(end); tvec1(end),yvec1(1); tvec1(1), yvec1(1)];
        
        
        figH = figure; subplot(2,1,1)
        imagesc(twin,xwin,xProj); colormap gray; axis image
        hold on
        plot(xbox1(:,1),xbox1(:,2),'y');
        plot(outTracks{i}(:,3),outTracks{i}(:,1),'r');
        hold off
        
        subplot(2,1,2)
        imagesc(twin,ywin,yProj); colormap gray; axis image
        hold on
        plot(ybox1(:,1),ybox1(:,2),'y');
        plot(outTracks{i}(:,3),outTracks{i}(:,2),'r');
        hold off
        
        answer = questdlg('Found another track within this box. What would you like to do with it?', 'Similar Track', 'Merge', 'Delete','Keep','Keep');
        
        if strcmp(answer,'Delete')
            removed = [removed; i];
        end
        if strcmp(answer,'Merge')
            merged = [merged; i];
        end
        close(figH);
    end
end
%             end

removed = unique(removed);
merged = unique(merged);

nRemoved = size(removed,1);



