function IM = refImage_ROIs(Results)

refIm = Results.Process.RefImage;

refIm_d = double(refIm);
refIm_d = (refIm_d - min(refIm_d(:)))./(max(refIm_d(:)) - min(refIm_d(:)));

IM_r = refIm_d;
IM_g = refIm_d;
IM_b = refIm_d;

ROI_IM = zeros(size(refIm_d));

for i = 1:size(Results.Process.ROIimage);
    curROI_IM = Results.Process.ROIimage{i,:};
    er = imerode(curROI_IM,strel('disk',1));
    
    curROI_bord = curROI_IM - er;
    ROI_IM(curROI_bord > 0) = 1;
end

IM_r(ROI_IM > 0) = 1;

IM = cat(3,IM_r,IM_g,IM_b);
