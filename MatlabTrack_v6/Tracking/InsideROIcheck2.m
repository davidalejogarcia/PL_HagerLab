function CentroidInAllRois = InsideROIcheck2(Centroid, ROIimage,varargin) 

if isempty(varargin)
    zero_flag = 1;
else
    zero_flag = varargin{1};
end
CentroidInRoi = cell(size(ROIimage));

Centroid = Centroid(Centroid(:,1) > 0,:);
NParticles = length(Centroid(:,1));

ROI_nest_mat = roiNest(ROIimage);



if size(Centroid,2) > 10
    xspec = 10;
    yspec = 11;
    ROIspec = 13;
    frameInd = 6;
else
    xspec = 1;
    yspec = 2;
    if size(Centroid,2) > 5
        ROIspec = 7;
        frameInd = 6;
    else
        ROIspec = 5;
        frameInd = 3;
    end
end
%Generate a background image where no ROIs are present
% BG = false(size(ROIimage{1}));
% for i = 1:length(ROIimage)
%     BG = BG + ROIimage{i};
%     
% end
%  BG = BG == 0;   
for ROIindx = 1:length(ROIimage)
    ROI_image_cur = ROIimage{ROIindx};
    isInside = zeros(NParticles,1);
    for ROIindx2 = 1:length(ROIimage)
        if ROI_nest_mat(ROIindx,ROIindx2) == 1
            ROI_image_cur = ROI_image_cur - ROIimage{ROIindx2};
            
        end
    end
%     ROI_image_cur(BG) = 0;
    
    for particleIx = 1:NParticles
        
        particlePos = Centroid(particleIx,:);
        if particlePos(xspec) > 0
            if ROI_image_cur(max(1,round(particlePos(yspec))),max(1,round(particlePos(xspec)))) == 1
                isInside(particleIx) = 1;
            end

    %         j = length(ROIpos(:,1));
    % 
    %         for i =1:1:length(ROIpos(:,1));
    % 
    %             if (ROIpos(i,2) < particlePos(2)) ~= (ROIpos(j,2) < particlePos (2)) ...
    %                     && ROIpos(i,1) + (particlePos(2)-ROIpos(i,2))/(ROIpos(j,2)-ROIpos(i,2))*...
    %                     (ROIpos(j,1)-ROIpos(i,1)) < particlePos(1);
    % 
    %                     isInside(particleIx) = not(isInside(particleIx));
    %             end
    %             j =i;
    %         end
        else
            if (ROIspec == 7 && size(particlePos,2) == 7) || (ROIspec == 13 && size(particlePos,2) == 13) || (ROIspec == 5 && size(particlePos,2) == 5) 
                isInside(particleIx) = 0;
            end
        end

    end
    Index = find(isInside ==1);
  
    
    CentroidInRoi{ROIindx} = Centroid(Index,:);
    if size(Centroid,2) > 5 && zero_flag == 1
        for i = 1: max(CentroidInRoi{ROIindx}(:,6))
            Index = find(CentroidInRoi{ROIindx}(:,6) == i);

            if isempty(Index);
                addvec = [0 0 0 0 0 i];
                if ROIspec == 13
                    addvec = [addvec zeros(1,7)];
                end
                if size(CentroidInRoi{ROIindx},2) == ROIspec
                    addvec(:,ROIspec) = ROIindx;
                end
                CentroidInRoi{ROIindx} = [CentroidInRoi{ROIindx}; addvec];
            end
        end
    end

     
     CentroidInRoi{ROIindx}(:,ROIspec) = ROIindx;
end
CentroidInAllRois = [];
for i = 1:length(ROIimage)
    CentroidInAllRois = [CentroidInAllRois; CentroidInRoi{i}];
end
%If we didn't find any particles at a particular time insert a dummy row
% if size(Centroid,2) > 5 && zero_flag == 1
%     for i = 1:max(CentroidInAllRois(:,6))
%         Index = find(CentroidInAllRoi(:,6) == i);
%         
%         if isempty(Index);
%             addvec = [0 0 0 0 0 i];
%             if ROIspec == 13
%                 addvec = [addvec zeros(1,7)];
%             end
%             
%             CentroidInRoi{ROIindx} = [CentroidInRoi{ROIindx}; addvec];
%         end
%     end
% end

if size(Centroid,2) > 5
    CentroidInAllRois = sortrows( CentroidInAllRois, [frameInd ROIspec]);
else
    CentroidInAllRois = sortrows(CentroidInAllRois, [4 frameInd]);
end
            
    

       
      
    