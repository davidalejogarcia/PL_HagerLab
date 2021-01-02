function  plotTwotrSegments(draw, trSegment1, trSegment2, axesHandle,varargin)

if isempty(varargin)
    plotSpatTime = 1;
     
else
    plotSpatTime = varargin{1};
    
end

switch draw
    case 1
        axes(axesHandle);
        delete(findobj(gca,'Color','g'));
        
        hold on;
        if plotSpatTime == 1
            plot(trSegment1(:,1),trSegment1(:,2),'r','LineWidth',1);
            plot(trSegment2(:,1),trSegment2(:,2),'y','LineWidth',1);
            hold off;
        elseif plotSpatTime == 3
            timevec1 = trSegment1(:,3);
            timevec2 = trSegment2(:,3);
            spatvec1 = trSegment1(:,1);
            spatvec2 = trSegment2(:,1);
            plot(timevec1,spatvec1,'r','LineWidth',1);
            plot(timevec2,spatvec2,'y','LineWidth',1);
            hold off;
        elseif plotSpatTime == 2
            timevec1 = trSegment1(:,3);
            timevec2 = trSegment2(:,3);
            spatvec1 = trSegment1(:,2);
            spatvec2 = trSegment2(:,2);
            plot(timevec1,spatvec1,'r','LineWidth',1);
            plot(timevec2,spatvec2,'y','LineWidth',1);
            hold off;
        end

    
    case 0
    delete(findobj(gca,'Color','r'));
    delete(findobj(gca,'Color','y'));
end