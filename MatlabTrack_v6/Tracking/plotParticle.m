function plotParticle(Particles,imageIx, isPSF,varargin)

% Plot particles using PSF fitting coordinates if isPSF == 1. Otherwise the
% Centroids coordinates are plot.

if isempty(varargin)
    ax_h = gca;
else
    ax_h = varargin{1};
end
delete(findobj(ax_h,'LineWidth',1));
origCO = get(groot,'DefaultAxesColorOrder');
origCO = repmat(origCO,[10,1]);
if isPSF
    ROIindx = 13;
else
    ROIindx = 7;
end
for i = 1:max(Particles(:,ROIindx))
    Particles_tmp = Particles(Particles(:,ROIindx) == i,:);
    Particles_tmp(Particles_tmp(:,1) == 0,:) = [];
    
    pIx = find(Particles_tmp(:,6) == imageIx);
    if ~isempty(pIx)
        hold on;
        if isPSF
            plot (Particles_tmp(pIx,10),Particles_tmp(pIx,11),'o', 'MarkerSize', 12,'MarkerEdgeColor',origCO(i,:));
        else
            plot (Particles_tmp(pIx,1),Particles_tmp(pIx,2),'o', 'MarkerSize', 12,'MarkerEdgeColor',origCO(i,:));
        end
        hold off;
    end
end
%Plot automatically added Particles in red
if ~isempty(Particles(Particles(:,ROIindx) == 0))
    Particles_added_tmp = Particles(Particles(:,ROIindx) == 0,:);
    pIx2 = find(Particles_added_tmp(:,6) == imageIx);
    hold on;
    if isPSF
        plot (Particles_added_tmp(pIx2,10),Particles_added_tmp(pIx2,11),'or', 'MarkerSize', 12);
    else
        plot (Particles_added_tmp(pIx2,1),Particles_added_tmp(pIx2,2),'or', 'MarkerSize', 12);
    end
    hold off;
end