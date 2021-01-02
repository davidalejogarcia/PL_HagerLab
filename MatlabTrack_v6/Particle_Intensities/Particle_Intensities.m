function Gmask_properties = Particle_Intensities(imStack,Tracks,varargin)

%Wrapper for gmask function that takes the positions from input Tracks, and
%cuts the appropriate part out of the correct image frame in ImStack, and
%returns the list of [x0,y0,I] for each particle.
%Intended for use with saved results generated from integratedTrack_GUI
%In addition to the image stack & track list, a structure of parameters can
%be specified. this structure can have any or all of the following fields:
% pxSize: Camera pixel size in microns in object space (default: 0.104)
% NA: Numerical Aperature of objective (default: 1.45)
% wavelength: emission wavelength of chromophore (default: 0.7)
% windowSize: size around spot used for fitting (default: 7)
%
%David Ball 11/2015

if ~isempty(varargin)
    if isfield(varargin{1},'pxSize');
        px_size = varargin{1}.pxSize;
    else
        px_size = 0.104;
    end
    if isfield(varargin{1},'NA');
        NA = varargin{1}.NA;
    else
        NA = 1.45;
    end    
    if isfield(varargin{1},'wavelength');
        Wl = varargin{1}.wavelength;
    else
        Wl = 0.7;
    end   
    if isfield(varargin{1},'windowSize')
        windowSize = varargin{1}.windowSize;
    else
        windowSize = 7;
    end
else
    
    px_size = 0.104;
    NA = 1.45;
    Wl = 0.7;
    windowSize = 7;
end
parameters.psfwidth = 2*0.61*Wl/NA/px_size;
% parameters.psfwidth = 10;
xdim = imStack(1).width;
ydim = imStack(1).height;


% Nframes = size(imStack,1);
Gmask_properties = [];

    
for j = 1:size(Tracks,1)
    curTrack = Tracks(j,:);
    im_array = imStack(curTrack(1,3)).data;
    
    
    xCenter = curTrack(1);
    yCenter = curTrack(2);
    
    xRound = round(xCenter);
    yRound = round(yCenter);
    if xRound > windowSize && yRound > windowSize && xRound < xdim - windowSize && yRound < ydim - windowSize
        imCut = im_array(yRound - windowSize:yRound + windowSize,xRound - windowSize: xRound + windowSize);
        x0 = xCenter - xRound + windowSize + 1;
        y0 = yCenter - yRound + windowSize + 1;
        results = gmask(imCut,parameters,x0,y0);
        results(1:2) = [xCenter,yCenter];
        Gmask_properties = [Gmask_properties; results];
    end
    
end


