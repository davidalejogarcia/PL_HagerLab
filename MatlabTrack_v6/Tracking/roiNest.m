function ROI_nest_matrix = roiNest(ROI_image)

%determines whether any of the ROIs are completely enclosed by other ROIs
% Input ROI_image should be a Nx1 cell array where each cell
% contains a binary image of an ROI
% Output is a NxN matrix of 0's and 1's where the (i,j) component tells if
% ROI{j} is completely interior to ROI{i}

nROIs = length(ROI_image);
ROI_nest_matrix = zeros(nROIs);

if nROIs > 1
    for i = 1:nROIs
        for j = 1:nROIs
            if j~=i
                ROI1 = ROI_image{i};
                ROI2 = ROI_image{j};
                ind1 = find(ROI1);
                ind2 = find(ROI2);
                
                if length(intersect(ind1,ind2)) == length(ind2);
                    ROI_nest_matrix(i,j) = 1;
                end
            end
            
        end
    end
end
for i = 1:size(ROI_nest_matrix,1)
    innerROIs = find(ROI_nest_matrix(i,:) == 1);
    if length(innerROIs) > 1
        for j = 1:length(innerROIs)
            for k = 1:length(innerROIs)
                if j ~= k && ROI_nest_matrix(innerROIs(j),innerROIs(k)) == 1
                    ROI_nest_matrix(i,innerROIs(k)) = 0;
                end
            end
        end
    end
end