function [Tracks,TrkPtsAdded] = ConvexStrict2FrameTracking(Particles,varargin)

%Particles should be a 3 or 4 column matrix, with the last column the
%time-points, and the first columns the spatial information

if ~isempty(varargin)
    parameters = varargin{1};
    if ~isfield(parameters,'Gv')
        parameters.Gv = 4;
    end
    if ~isfield(parameters,'Gd')
        parameters.Gd = 4;
    end
    if ~isfield(parameters,'shTr')
        parameters.shTr = 6;
    end
    if ~isfield(parameters,'gaps')
        parameters.gaps = 3;
    end
    
else
    parameters.Gv = 8;
    parameters.Gd = 8;
    parameters.shTr = 3;
    parameters.gaps = 3;
end

%determine the number of spatial dimensions
dim = size(Particles,2);
sp_dim = dim - 1;

tpoints = unique(Particles(:,dim));

%Get the particles in the first existing timepoint & use these as a
%starting point for the Tracks that will be output
T = Particles(Particles(:,dim) == tpoints(1),1:sp_dim);
Tracks = [T,tpoints(1)*ones(size(T,1),1),(1:size(T,1))'];

%add a place holder for the velocity,timepoints in velocity calculation,
%cost, the last timepoint the particle was found, and a flag on if the
%particle was lost
T = [T, NaN(size(T)),zeros(size(T,1),2),tpoints(1)*ones(size(T,1),1),ones(size(T,1),1)];
maxInd = size(T,1);
for nIter = 2:size(tpoints,1)

    Pi = Particles(Particles(:,dim) == tpoints(nIter),:);

    costMat = costCalc(T,Pi,parameters.Gd,parameters.Gv,parameters.gaps);

    [assign,totalCost] = assignmentoptimal(costMat);
   
    
    
    %Reassign values in T based on current track points
    for i = 1:size(assign,1)
        if assign(i) ~= 0
            T(assign(i),sp_dim+1:2*sp_dim) = (Pi(i,1:sp_dim) - T(assign(i),1:sp_dim))/(tpoints(nIter)-T(assign(i),2*sp_dim+3));
            T(assign(i),1:sp_dim) = Pi(i,1:sp_dim);
            T(assign(i),2*sp_dim+1) = T(assign(i),2*sp_dim+1) + 1;
            T(assign(i),2*sp_dim+2) = T(assign(i),2*sp_dim+2) + costMat(i,assign(i));
            T(assign(i),2*sp_dim+3) = tpoints(nIter);
        end
    end
    
    %Remove unsused tracks from further analysis
    unusedTracks = ones(size(T,1),1);
    unusedTracks(assign(assign > 0)) = 0;
    
    T(unusedTracks == 1,end) = 0;
    
    %add newly identified points to the tracks to check
    Tnew = [Pi(assign == 0,1:sp_dim),NaN(size(Pi(assign == 0,1:sp_dim))),zeros(size(Pi(assign == 0,:),1),2),tpoints(nIter)*ones(size(Pi(assign == 0,:),1),1),ones(size(Pi(assign == 0,:),1),1)];
    T = [T;Tnew];
    
    %Add new track indices to particles that couldn't be assigned
    for i = 1:size(assign,1)
        if assign(i) == 0
            maxInd = maxInd + 1;
            assign(i) = maxInd;
        end
    end
    
    
    TracksNew = [Pi,assign];
    Tracks = [Tracks;TracksNew];
end

Tracks = sortrows(Tracks,[4,3]);


%  Deleting 'zeros' lines 
% (correspoding to frames where no track has beenfound


Tracks(Tracks(:,1)==0,:)= [];


% Fill Gaps (blinking) in tracks
TrkPtsAdded = cell(max(Tracks(:,4)),1);
for i = 1: max(Tracks(:,4))
    
    tempTrack = Tracks(Tracks(:,4) == i,:);
    diff = tempTrack(2:end,3)-tempTrack(1:end-1,3);
    lIx = find (diff ~=1);
    if  ~isempty(lIx);
        for iGap = 1:length(lIx)
            
            for iFrame = 1:diff(lIx(iGap))-1
                catTrack(:,1) =  tempTrack(lIx(iGap),1) + iFrame*(tempTrack(lIx(iGap)+1,1) - tempTrack(lIx(iGap),1))/diff(lIx(iGap));
                catTrack(:,2) = tempTrack(lIx(iGap),2)+ iFrame*(tempTrack(lIx(iGap)+1,2) - tempTrack(lIx(iGap),2))/diff(lIx(iGap));
                catTrack(:,3) = tempTrack(lIx(iGap),3) + iFrame;%round((tempTrack(lIx+1,3) + tempTrack(lIx,3))/diff(lIx));
                catTrack(:,4) = tempTrack(lIx(iGap),4);         
                Tracks = [Tracks;catTrack];
                %keep a note of which particle points were added
                TrkPtsAdded{i} = [TrkPtsAdded{i}; tempTrack(lIx(iGap),3) + iFrame];
            end
        end
        clear catTrack;
    end
    
%     TrkPtsAdded{i} = (tempTrack(lIx+1,3) + tempTrack(lIx,3))/2;
end



% Re-sort track with the filled gaps
Tracks=sortrows(Tracks, [4 3]);

%remove tracks that are shorter than the specified shortest track
ind = 1;
Tracks2 = [];

for i = 1:max(Tracks(:,4))
    tmp = Tracks(Tracks(:,4) == i,:);
    if size(tmp,1) >= parameters.shTr
        tmp(:,4) = ind;
        ind = ind+1;
        Tracks2 = [Tracks2;tmp];
    end
    
end
Tracks = Tracks2;


function costMat = costCalc(T1,Pi,Gd,Gv,gaps)
%convex cost function taking into account displacement & velocity

%Determine the number of spatial dimensions
ndimT1 = (size(T1,2) - 4)/2;

%Calculate the weights for the convex combination
vel_vecs = T1(:,ndimT1+1:ndimT1+ndimT1);
vel_norms = sqrt(sum((vel_vecs).^2,2));

mu_D = exp(-vel_norms);
mu_D(isnan(mu_D)) = 1;

mu_V = exp(-1*((vel_norms - Gv/2).^2)/2);
mu_V(isnan(mu_V)) = 0;


%construct the cost matrix

for i = 1:size(Pi,1)
    for j = 1:size(T1,1)
        costMat(i,j) = sqrt(sum((Pi(i,1:ndimT1) - T1(j,1:ndimT1)).^2));
        if costMat(i,j) > Gd || (Pi(i,ndimT1+1) - T1(j,2*ndimT1+3)) > gaps || T1(j,2*ndimT1+4) == 0
            costMat(i,j) = Inf;
        end
        costMat(i,j) = mu_D(j,:)*costMat(i,j);
        
        vel_mot = T1(j,1:ndimT1) + T1(j,ndimT1+1:2*ndimT1)*(Pi(i,ndimT1+1) - T1(j,2*ndimT1+3));
        
        costMatV(i,j) = mu_V(j,:)*sqrt(sum((Pi(i,1:ndimT1) - vel_mot).^2));
        if ~isnan(costMatV(i,j))
            costMat(i,j) = costMat(i,j) + costMatV(i,j);
        end
    end
end


        
