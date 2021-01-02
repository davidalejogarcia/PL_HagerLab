function plotSingleTrack (Track, imageIx, axesHandle,varargin)

if isempty(varargin)
    plotSpatTime = 1;
     
else
    plotSpatTime = varargin{1};
    
end
axes(axesHandle);
delete(findobj(gca,'Color','g'));
        
        hold on;
        if plotSpatTime == 1
            plotIx = find(Track(:,3) <= imageIx);
            plot(Track(plotIx,1),Track(plotIx,2),'g','LineWidth',1);
            hold off;
        elseif plotSpatTime == 2
            timevec = Track(:,3);
            spatvec = Track(:,1);
            plot(timevec,spatvec,'g','LineWidth',1);
            hold off;
        elseif plotSpatTime == 3
            timevec = Track(:,3);
            spatvec = Track(:,2);
            plot(timevec,spatvec,'g','LineWidth',1);
            hold off;
        end