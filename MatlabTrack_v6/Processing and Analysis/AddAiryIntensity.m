function OutTracks = AddAiryIntensity(InTracks,pixelSize, Stack)

Wl = 0.700; % Wavelength in microns
NA = 1.45;  % Numerical aperture of the objective

radiusAiry = 2*0.61*Wl/NA/pixelSize;
radiusAiry = 3;
NPeaks = length(InTracks(:,1));
[height,width] = size(Stack(1).data);

[meshX, meshY] = meshgrid(1:width, 1:height);




for i = 1:NPeaks
    Peak = InTracks(i,:);
    Frame = Peak(3);
    if i == 2359
        aia = 1;
    end
    CheckMatrix = (Peak(1) - meshX).^2 + (Peak(2)-meshY).^2;
    idx = find (CheckMatrix<= radiusAiry^2);
    Npixel = length(idx);
    
    Intensity = sum(Stack(Frame).data(idx))/Npixel;
    
    idx = find (CheckMatrix > radiusAiry^2 & CheckMatrix <= (radiusAiry + 1)^2);
    Npixel = length(idx);
    bkg = sum(Stack(Frame).data(idx))/Npixel;
    
    OutTracks(i,:) = [Peak, Intensity, bkg];
   
end  

    
