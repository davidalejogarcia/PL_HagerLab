function distance = distanceFromPol(measCentr,refCentr)

% Find distance between measurement and reference spots
distance = [];
for i = 1:length(measCentr(:,1));
    xTemp = measCentr(i,1);
    yTemp = measCentr(i,2);
    distance(i) = min(sqrt((refCentr(:,1)-xTemp).^2 +(refCentr(:,2)-yTemp).^2));
    
end