function CentroidInRoi = InsideROIcheck(Centroid, ROIpos) 



NParticles = length(Centroid(:,1));
isInside = zeros(NParticles,1);

for particleIx = 1:NParticles
    
    particlePos = Centroid(particleIx,:);
    
   
    j = length(ROIpos(:,1));
    
    for i =1:1:length(ROIpos(:,1));
    
        if (ROIpos(i,2) < particlePos(2)) ~= (ROIpos(j,2) < particlePos (2)) ...
                && ROIpos(i,1) + (particlePos(2)-ROIpos(i,2))/(ROIpos(j,2)-ROIpos(i,2))*...
                (ROIpos(j,1)-ROIpos(i,1)) < particlePos(1);
            
                isInside(particleIx) = not(isInside(particleIx));
        end
    j =i;
    end

end
    Index = find(isInside ==1);
  
    
    CentroidInRoi = Centroid(Index,:);
    
    for i = 1: max(CentroidInRoi(:,6))
        Index = find(CentroidInRoi(:,6) == i);
        
        if isempty(Index);
            CentroidInRoi = [CentroidInRoi; 0 0 0 0 0 i];
        end
    end
    
     CentroidInRoi=sortrows( CentroidInRoi, 6);
    
            
    

       
      
    