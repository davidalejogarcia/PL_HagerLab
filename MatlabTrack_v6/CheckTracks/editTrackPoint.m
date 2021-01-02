function newPos = editTrackPoint(axisHandle)
   
    
% Click on the new position of the particle

h_point = impoint(axisHandle);
newPos = getPosition(h_point);




end
