function plotROI_t(ROI,dim,tvec,varargin)

if isempty(varargin)
    ax_h = gca;
    selROI = 0;
elseif length(varargin) == 1
    ax_h = varargin{1};
    selROI = 0;
else
    ax_h = varargin{1};
    selROI = varargin{2};
end


delete(findobj(ax_h,'LineWidth',1.1));
origCO = get(groot,'DefaultAxesColorOrder');
origCO = repmat(origCO,[100,1]);

axes(ax_h);
for i = 1:length(ROI)
    ROI2plot = ROI{i,:}(:,dim);
    vec1 = [min(ROI2plot),min(ROI2plot)];
    vec2 = [max(ROI2plot),max(ROI2plot)];
    
    hold on;
    if i == selROI
        plot(tvec,vec1,'LineWidth',1.1,'Color',[1 1 1]);
        plot(tvec,vec2,'LineWidth',1.1,'Color',[1 1 1]);
    else
        plot(tvec,vec1,'LineWidth',1.1,'Color',origCO(i,:));
        plot(tvec,vec2,'LineWidth',1.1,'Color',origCO(i,:));
    end
    hold off;
end