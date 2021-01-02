function [Tracks,TrkPtsAdded] = ngaTracking(Particles,varargin)

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
    if ~isfield(parameters,'Twin')
        parameters.Twin = 3;
    end
    if ~isfield(parameters,'shTr')
        parameters.shTr = 2;
    end
    if ~isfield(parameters,'gaps')
        parameters.gaps = 3;
    end
else
    parameters.Gv = 8;
    parameters.Gd = 8;
    parameters.Twin = 3;
    parameters.shTr = 2;
    parameters.gaps = 3;
end

%determine the number of spatial dimensions
dim = size(Particles,2);
sp_dim = dim - 1;

tpoints = (min(Particles(:,dim)):max(Particles(:,dim)))';
assign_all = [];
%Separate the particles in each frame
PartFrames = cell(size(tpoints,1),1);
for i = 1:size(tpoints,1)
    PartFrames{i,:} = Particles(Particles(:,dim) == tpoints(i),:);
end

%Add track numbers to the first timepoint
PartFrames{1,:}(:,dim+1) = (1:size(PartFrames{1,:},1))';
maxTrInd = size(PartFrames{1,:},1);

%Get the particles in the first existing timepoint & use these as a
%starting point for the Tracks that will be output
costMatOrig = cell(parameters.Twin - 1,1);
assignOrig = cell(parameters.Twin - 1,1);
T = cell(size(tpoints,1),1);

T{1,:} = PartFrames{1,:}(:,1:sp_dim);
%column order: 1-sp_dim: spatial coordinates of last known location
%              sp_dim+1:2*sp_dim: velocity components
%              2*sp_dim+1: number of points used to calcualate an avg.
%              velcity (not presently used)
%              2*sp_dim+2: cost value of last connection
%              2*sp_dim+3: last frame where the track was found
%              2*sp_dim+4: flag (0 or 1) if track should be considered for
%              extension
%            
T{1,:} = [T{1,:}, NaN(size(T{1,:})),zeros(size(T{1,:},1),2),tpoints(1)*ones(size(T{1,:},1),1),ones(size(T{1,:},1),1)];

for i = 2:parameters.Twin - 1
    T{i,:} = T{i-1,:};
    
    Pi = PartFrames{i,:};
    costMatOrig{i,:} = costCalc(T{i,:},Pi,parameters.Gd,parameters.Gv,parameters.gaps);
    [assign,~] = assignmentoptimal(costMatOrig{i,:});
    
    costVec = zeros(size(assign));
    for ind = 1:size(assign,1)
        if assign(ind,:) > 0
                costVec(ind) = costMatOrig{i,:}(ind,assign(ind,:));
            else
                costVec(ind) = Inf;
            end
%         costVec(ind) = costMatOrig{i,:}(ind,assign(ind));
    end
    assignOrig{i,:} = [assign,costVec];
    %Reassign values in T based on current track points
    for j = 1:size(assign,1)
        if assign(j) ~= 0
            T{i,:}(assign(j),sp_dim+1:2*sp_dim) = (Pi(j,1:sp_dim) - T{i,:}(assign(j),1:sp_dim))/(tpoints(i+1)-T{i,:}(assign(j),end));
            T{i,:}(assign(j),1:sp_dim) = Pi(j,1:sp_dim);
            T{i,:}(assign(j),2*sp_dim+1) = T{i,:}(assign(j),2*sp_dim+1) + 1;
            T{i,:}(assign(j),2*sp_dim+2) = T{i,:}(assign(j),2*sp_dim+2) + costMatOrig{i,:}(j,assign(j));
            T{i,:}(assign(j),2*sp_dim+3) = tpoints(i+1);
        end
    end
    %add newly identified points to the tracks to check
    Tnew = [Pi(assign == 0,1:sp_dim),NaN(size(Pi(assign == 0,1:sp_dim))),zeros(size(Pi(assign == 0,:),1),2),tpoints(i)*ones(size(Pi(assign == 0,:),1),1),ones(size(Pi(assign == 0,:),1),1)];
    T{i,:} = [T{i,:};Tnew];

    
    %Add new track indices to particles that couldn't be assigned
    for j = 1:size(assign,1)
        if assign(j) == 0
            maxTrInd = maxTrInd + 1;
            assign(j) = maxTrInd;
        end
    end
    PartFrames{i,:}(:,dim+1) = assign;
    
    
    
end

for nIter = parameters.Twin:size(tpoints,1)
    T{nIter,:} = T{nIter - 1,:};
    if nIter == 9
        na = 10;
    end
    Pi = Particles(Particles(:,dim) == tpoints(nIter),:);
    
    %compare frames: F(i-Twin+1) and F(i), then F(i-Twin+2) and F(i), ...
    %F(i-1) and F(i)
    costMat = cell(parameters.Twin-1,1);
    assignCost = zeros(size(Pi,1),2,parameters.Twin-1);
    assign = zeros(size(Pi,1),parameters.Twin-1);
    for iWin = nIter - parameters.Twin + 1: nIter - 1
        Tcur = T{iWin,:};
        
        costMat{iWin-nIter+parameters.Twin,:} = costCalc(Tcur,Pi,parameters.Gd,parameters.Gv,parameters.gaps);

        [assign(:,iWin-nIter+parameters.Twin),~] = assignmentoptimal(costMat{iWin-nIter+parameters.Twin,:});
        costVec = zeros(size(assign(:,iWin-nIter+parameters.Twin)));
        
        for ind = 1:size(assign,1)
            if assign(ind,iWin-nIter+parameters.Twin) > 0
                costVec(ind) = costMat{iWin-nIter+parameters.Twin,:}(ind,assign(ind,iWin-nIter+parameters.Twin));
            else
                costVec(ind) = Inf;
            end
        end
        assignCost(:,:,iWin-nIter+parameters.Twin) = [assign(:,iWin-nIter+parameters.Twin),costVec];
    end
    
    
    for i = 1:size(assign,1)
        lastIndCur = unique(assign(i,assign(i,:) > 0));
        if isempty(lastIndCur)
            maxTrInd = maxTrInd + 1;
            PartFrames{nIter,:}(i,dim+1) = maxTrInd;
            Tnew = [Pi(i,1:sp_dim),NaN(size(Pi(i,1:sp_dim))),zeros(size(Pi(i,:),1),2),tpoints(nIter)*ones(size(Pi(i,:),1),1),ones(size(Pi(i,:),1),1)];
            T{nIter,:} = [T{nIter,:};Tnew];
        elseif length(lastIndCur) == 1
            PartFrames{nIter,:}(i,dim+1) = lastIndCur;
            T{nIter,:}(lastIndCur,sp_dim+1:2*sp_dim) = (Pi(i,1:sp_dim) - T{nIter,:}(lastIndCur,1:sp_dim))/(tpoints(nIter)-T{nIter,:}(lastIndCur,end));
            T{nIter,:}(lastIndCur,1:sp_dim) = Pi(i,1:sp_dim);
            T{nIter,:}(lastIndCur,2*sp_dim+1) = T{nIter,:}(lastIndCur,2*sp_dim+1) + 1;
            T{nIter,:}(lastIndCur,2*sp_dim+2) = T{nIter,:}(lastIndCur,2*sp_dim+2) + costMat{parameters.Twin-1,:}(i,lastIndCur);
            T{nIter,:}(lastIndCur,2*sp_dim+3) = tpoints(nIter);
        elseif length(lastIndCur) > 1
            
            TotalCost = zeros(length(lastIndCur),1);
            for j = 1:length(lastIndCur)
                assignInd = find(assign(i,:) == lastIndCur(j));
%                 if assign(assignInd) ~= 0
                    IndCost = assignCost(i,2,assignInd);
                    IndCost = IndCost(:);
                    IndCostRecip = 1./IndCost;
                
                    TotalCost(j) = 1./(sum(IndCostRecip));
%                 else
%                     TotalCost(j) = Inf;
%                 end
            end
            [~,MinInd] = min(TotalCost);
            lastIndCur = lastIndCur(MinInd);
            WinInd = find(assign(i,:) == lastIndCur,1,'last');
            PartFrames{nIter,:}(i,dim+1) = lastIndCur;
            if nIter == 3 && i == 3
                noa = 1;
            end
            
            T{nIter,:}(lastIndCur,sp_dim+1:2*sp_dim) = (Pi(i,1:sp_dim) - T{nIter,:}(lastIndCur,1:sp_dim))/(tpoints(nIter)-T{nIter,:}(lastIndCur,end));
            T{nIter,:}(lastIndCur,1:sp_dim) = Pi(i,1:sp_dim);
            T{nIter,:}(lastIndCur,2*sp_dim+1) = T{nIter,:}(lastIndCur,2*sp_dim+1) + 1;
            T{nIter,:}(lastIndCur,2*sp_dim+2) = T{nIter,:}(lastIndCur,2*sp_dim+2) + costMat{WinInd,:}(i,lastIndCur);
                        
            
            
        end
        rem_tpoint = nIter - parameters.gaps - 1;
        T{nIter,:}(T{nIter,:}(:,2*sp_dim+3) == rem_tpoint,2*sp_dim+4) = 0;
        
    end
        
            
%     %Reassign values in T based on current track points
%     for i = 1:size(assign,2)-1
%         for j = 1:size(assign,1)
%             
%             if assign(j,size(assign,2)-i+1) ~= 0
%                 T(assign(i),sp_dim+1:2*sp_dim) = (Pi(i,1:sp_dim) - T(assign(i),1:sp_dim))/(tpoints(nIter)-T(assign(i),end));
%                 T(assign(i),1:sp_dim) = Pi(i,1:sp_dim);
%                 T(assign(i),2*sp_dim+1) = T(assign(i),2*sp_dim+1) + 1;
%                 T(assign(i),2*sp_dim+2) = T(assign(i),2*sp_dim+2) + costMat(i,assign(i));
%                 T(assign(i),2*sp_dim+3) = tpoints(nIter);
%             end
%         end
%     end
%     %add newly identified points to the tracks to check
%     Tnew = [Pi(assign == 0,1:sp_dim),NaN(size(Pi(assign == 0,1:sp_dim))),zeros(size(Pi(assign == 0,:),1),2),tpoints(nIter)*ones(size(Pi(assign == 0,:),1),1)];
%     T = [T;Tnew];
%     %Add new track indices to particles that couldn't be assigned
%     for i = 1:size(assign,1)
%         if assign(i) == 0
%             maxInd = maxInd + 1;
%             assign(i) = maxInd;
%         end
%     end
%     
%     
%     TracksNew = [Pi,assign];
%     Tracks = [Tracks;TracksNew];
end
Tracks = [];
for i = 1:size(PartFrames)
    Tracks = [Tracks;PartFrames{i,:}];
    
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