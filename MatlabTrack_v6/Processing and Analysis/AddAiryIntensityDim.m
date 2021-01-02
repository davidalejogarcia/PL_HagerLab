function OutIntensity = AddAiryIntensityDim(InTracks,pixelSize, Stack)

Wl = 0.600; % Wavelength in microns
NA = 1.45;  % Numerical aperture of the objective

radiusAiry = 0.61*Wl/NA/pixelSize;
NTracks = length(InTracks(:,1));
[height,width] = size(Stack(1).data);

[meshX, meshY] = meshgrid(1:width, 1:height);


OutIntensity = [];

for i = 1:NTracks
    
    curTrack = InTracks(i,:);
    OutIntensity_cur = zeros(curTrack(1,2),3);
    for j = curTrack(1,5):curTrack(1,2)
        
        OutIntensity_cur(:,1) = curTrack(1,1);
        if curTrack(1,4) == 1
            xwin = curTrack(1,6):curTrack(1,6)+curTrack(1,3);
            ywin = curTrack(1,7):curTrack(1,8);
            xbg = [curTrack(1,6)-1:curTrack(1,6), curTrack(1,6)+curTrack(1,3):curTrack(1,6)+curTrack(1,3)+1];
            ybg = [curTrack(1,7)-1:curTrack(1,7),curTrack(1,8):curTrack(1,8)+1];
        else
            ywin = curTrack(1,6):curTrack(1,6)+curTrack(1,3);
            xwin = curTrack(1,7):curTrack(1,8);
            ybg = [curTrack(1,6)-1:curTrack(1,6), curTrack(1,6)+curTrack(1,3):curTrack(1,6)+curTrack(1,3)+1];
            xbg = [curTrack(1,7)-1:curTrack(1,7),curTrack(1,8):curTrack(1,8)+1];
        end
        PartIm = Stack(j).data(ywin,xwin);
        Npixel = numel(PartIm);
        
        
        OutIntensity_cur(j,2) = sum(PartIm(:))/Npixel;
        
        [bg_subx, bg_suby] = meshgrid(xbg,ybg);
        bg_subx = bg_subx(:);
        bg_suby = bg_suby(:);
        
        bg_ind = sub2ind(bg_suby,bg_subx);
        
        OutIntensity_cur(j,3) = sum(Stack(j).data(bg_ind))/numel(bg_ind);
    end
   OutIntensity = [OutIntensity;OutIntensity_cur];
end  

    
