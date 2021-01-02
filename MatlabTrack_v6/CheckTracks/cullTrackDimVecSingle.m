function [outTrackDim,removed] = cullTrackDimVecSingle(inTrackDim,curIx,wholeImage)

outTrackDim = inTrackDim;

%create a 3D image from the input stack
% for i = 1:length(imStack)
%     wholeImage(:,:,i) = imStack(i).data;
% end

    removed = [];
    for i = 1:size(inTrackDim,1)
        if isempty(find(removed == i,1,'first'))&& max(inTrackDim(curIx,:)) > 0 && i ~= curIx
            tvec1 = (inTrackDim(i,5):inTrackDim(i,5)+inTrackDim(i,2))';
            if inTrackDim(i,4) == 1
                xvec1 = (inTrackDim(i,6):inTrackDim(i,6)+inTrackDim(i,3))';
                yvec1 = (inTrackDim(i,7):inTrackDim(i,8))';
            else
                yvec1 = (inTrackDim(i,6):inTrackDim(i,6)+inTrackDim(i,3))';
                xvec1 = (inTrackDim(i,7):inTrackDim(i,8))';
            end
%             for j = i+1:size(inTrackDim,1) 
                if isempty(find(removed == curIx,1,'first')) && isempty(find(removed == i,1,'first'))
                    tvec2 = (inTrackDim(curIx,5):inTrackDim(curIx,5)+inTrackDim(curIx,2))';
                    if inTrackDim(curIx,4) == 1
                        xvec2 = (inTrackDim(curIx,6):inTrackDim(curIx,6)+inTrackDim(curIx,3))';
                        yvec2 = (inTrackDim(curIx,7):inTrackDim(curIx,8))';
                    else
                        yvec2 = (inTrackDim(curIx,6):inTrackDim(curIx,6)+inTrackDim(curIx,3))';
                        xvec2 = (inTrackDim(curIx,7):inTrackDim(curIx,8))';
                    end
                    
                    if (size(intersect(tvec1,tvec2),1)/min(size(tvec1,1),size(tvec2,1)) > 0.5) && (size(intersect(xvec1,xvec2),1)/min(size(xvec1,1),size(xvec2,1)) > 0.5) && (size(intersect(yvec1,yvec2),1)/min(size(yvec1,1),size(yvec2,1)) > 0.5)
                        twin = (min(tvec1(1),tvec2(1)):max(tvec1(end),tvec2(end)))';
                        curImage = wholeImage(:,:,twin);
                        xwin = (min(xvec1(1),xvec2(1)):max(xvec1(end),xvec2(end)))';
                        ywin = (min(yvec1(1),yvec2(1)):max(yvec1(end),yvec2(end)))';
                        xProj = max(permute(curImage(ywin,xwin,:),[2,3,1]),[],3);
                        yProj = max(permute(curImage(ywin,xwin,:),[1,3,2]),[],3);
                        
                        xbox1 = [tvec1(1),xvec1(1); tvec1(1),xvec1(end); tvec1(end),xvec1(end); tvec1(end),xvec1(1); tvec1(1), xvec1(1)];
                        ybox1 = [tvec1(1),yvec1(1); tvec1(1),yvec1(end); tvec1(end),yvec1(end); tvec1(end),yvec1(1); tvec1(1), yvec1(1)];
                        xbox2 = [tvec2(1),xvec2(1); tvec2(1),xvec2(end); tvec2(end),xvec2(end); tvec2(end),xvec2(1); tvec2(1), xvec2(1)];
                        ybox2 = [tvec2(1),yvec2(1); tvec2(1),yvec2(end); tvec2(end),yvec2(end); tvec2(end),yvec2(1); tvec2(1), yvec2(1)];
                        
                        figH = figure; subplot(2,1,1)
                        imagesc(twin,xwin,xProj); colormap gray
                        hold on
                        plot(xbox1(:,1),xbox1(:,2),'r');
                        plot(xbox2(:,1),xbox2(:,2),'y');
                        hold off
                        
                        subplot(2,1,2)
                        imagesc(twin,ywin,yProj); colormap gray
                        hold on
                        plot(ybox1(:,1),ybox1(:,2),'r');
                        plot(ybox2(:,1),ybox2(:,2),'y');
                        hold off
                        
                        answer = questdlg('Which Track do you want to keep?', 'Merge Tracks','Red','Yellow','Both','Both');
                        
                        switch answer
                            case 'Red'
                                removed = [removed; curIx];
                            case 'Yellow'
                                removed = [removed; i];
                        end
                        close(figH);
                    end
                end
%             end
        end
    end
    removed = unique(removed);
    outTrackDim(removed,:) = [];
    nRemoved = size(removed,1);
    inTrackDim = outTrackDim;


