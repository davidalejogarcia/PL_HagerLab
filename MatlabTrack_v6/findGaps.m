function [gapDist, gaps2close] = findGaps(varargin)
files = uigetfile('*.mat','Select the images to open','','MultiSelect','on');
savefold = uigetdir(pwd,'Select a save location');


%change this value to adjust the percentile
if isempty(varargin)
    thresh = 0.75;
else
    thresh = varargin{1};
end


gapDist = [];
if ~iscell(files)
    files = {files};
end

if ~isempty(files) 
    for i = 1:length(files)
        load(files{i});
        particles = Results.Tracking.Particles(Results.Tracking.Particles(:,6) <= 100,:);

        partInd = 1;
        usedInd = [];
        while partInd <= size(particles,1)
            if i == 31
                bal = 10;
            end
    %         partInd
            newPart = particles(partInd,10:11);
            if newPart(1) > 0
                dX = sqrt((particles(:,10) - newPart(1)).^2 + (particles(:,11) - newPart(2)).^2);
                ind1 = find(dX <= 1);
                usedInd = [usedInd; ind1];
                t = particles(ind1,6);
                if length(t) > 1
                    dT = diff(t);
                    gaps = dT(dT > 1);
                    gapDist = [gapDist; (gaps - 1)];
                end
            end
            partInd = partInd + 1;
            while ~isempty(find(usedInd == partInd,1,'first')) &&  partInd <= length(particles)
                partInd = partInd + 1;
            end
        end
    end

    [n, x] = hist(gapDist,1:max(gapDist));
    N = n./sum(n);
    C = cumsum(N);
    ind = find(C >= thresh,1,'first');

    gaps2close = x(ind);


    save([savefold, filesep,'gapDist.mat'],'gapDist','n','x','gaps2close');
end
    