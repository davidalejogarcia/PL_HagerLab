function [fractions,resTimes] = PieChartsFromResTimeAnalysis(Results,varargin)

    
if Results.FitPar.Surv_PB(1,8) < 0.05
    fractions(:,1) = [1.0-Results.FitPar.Surv_PB(1,6); Results.FitPar.Surv_PB(1,6)*Results.FitPar.Surv_PB(1,5); Results.FitPar.Surv_PB(1,6)*(1.0-Results.FitPar.Surv_PB(1,5))];

    fractions(:,2) = [Results.FitPar.Surv_PB(2,6); ...
        fractions(2,1)*sqrt((Results.FitPar.Surv_PB(2,5)/Results.FitPar.Surv_PB(1,5)).^2 +...
        (Results.FitPar.Surv_PB(2,6)/Results.FitPar.Surv_PB(1,6)).^2); ...
        fractions(2,1)*sqrt((Results.FitPar.Surv_PB(2,5)/Results.FitPar.Surv_PB(1,5)).^2 +...
        (Results.FitPar.Surv_PB(2,6)/Results.FitPar.Surv_PB(1,6)).^2)];

    resTimes(:,1) = [1/Results.FitPar.Surv_PB(1,3); 1/Results.FitPar.Surv_PB(1,4)];
    resTimes(:,2) = [(resTimes(1,1)/Results.FitPar.Surv_PB(1,3))*Results.FitPar.Surv_PB(2,3); ...
        (resTimes(2,1)/Results.FitPar.Surv_PB(1,4))*Results.FitPar.Surv_PB(2,4)];
%     resTimes(:,2) = [(1/Results.FitPar.Surv_PB(1,3)) - (1/(Results.FitPar.Surv_PB(1,3)+Results.FitPar.Surv_PB(2,3))) ; ...
%         (1/Results.FitPar.Surv_PB(1,4)) - (1/(Results.FitPar.Surv_PB(1,4)+Results.FitPar.Surv_PB(2,4)))];

    figure; pie(fractions(:,1),{'','',''})

    set(gcf,'Units','centimeters')
    set(gcf,'paperUnits','centimeters')

    pos = get(gcf,'Position');

    pos(3) = 17.8/2;
    pos(4) = 23.0/4;

    set(gcf,'Position',pos);
    set(gcf,'paperPosition',pos);
    pos_ax2 = [0.025, 0.015, 0.62,0.9];
    set(gca,'Position',pos_ax2)

    %May need to ajust the x/y positions of each text 
    %For real Data:
    text(1.1,-0.25,[num2str(round(fractions(2,1)*100,1)), '\pm', num2str(round(fractions(2,2)*100,1)),'%'],'FontSize',10)
    text(0.9,-0.5,['\tau = ',num2str(round(resTimes(1,1),3)), '\pm', num2str(round(resTimes(1,2),3)),' s'],'FontSize',10)
    
    text(0.25,1.35,[num2str(round(fractions(3,1)*100,1)), '\pm', num2str(round(fractions(3,2)*100,1)),'%'],'FontSize',10)
    text(0.05, 1.1,['\tau = ',num2str(round(resTimes(2,1),3)), '\pm', num2str(round(resTimes(2,2),3)),' s'],'FontSize',10)
    
    text(-1.2,1.15,[num2str(round(fractions(1,1)*100,1)), '\pm', num2str(round(fractions(1,2)*100,1)),'%'],'FontSize',10)
%     text(-0.85,0.9,[num2str(round(fractions(1,1)*100,1)), '\pm', num2str(round(fractions(1,2)*100,1)),'%'],'FontSize',10)

    %For simulated Data:
%     text(-0.35,-0.25,[num2str(round(fractions(1,1)*100,1)), '\pm', num2str(round(fractions(1,2)*100,1)),'%'],'FontSize',10,'Color','w')
% 
% 
%     text(0.9,0.85,[num2str(round(fractions(2,1)*100,1)), '\pm', num2str(round(fractions(2,2)*100,1)),'%'],'FontSize',10)
%     text(0.9, 0.6,['\tau = ',num2str(round(resTimes(1,1),3)), '\pm', num2str(round(resTimes(1,2),3)),' s'],'FontSize',10)
% 
%     % text(-0.6,1.15,[num2str(round(fractions(1,1)*100,1)), '\pm', num2str(round(fractions(1,2)*100,1)),'%'],'FontSize',10)
%     text(-0.35,1.35,[num2str(round(fractions(3,1)*100,1)), '\pm', num2str(round(fractions(3,2)*100,1)),'%'],'FontSize',10)
%     text(-0.55, 1.1,['\tau = ',num2str(round(resTimes(2,1),3)), '\pm', num2str(round(resTimes(2,2),3)),' s'],'FontSize',10)
else
    fractions(:,1) = [1.0-Results.FitPar.Surv_PB(1,2); Results.FitPar.Surv_PB(1,2)];

    fractions(:,2) = [Results.FitPar.Surv_PB(2,2); Results.FitPar.Surv_PB(2,2)];

    resTimes(:,1) = 1/Results.FitPar.Surv_PB(1,1);
    resTimes(:,2) = (resTimes(1,1)/Results.FitPar.Surv_PB(1,1))*Results.FitPar.Surv_PB(2,1);
%     resTimes(:,2) = (1/Results.FitPar.Surv_PB(1,1)) - (1/(Results.FitPar.Surv_PB(1,1)+Results.FitPar.Surv_PB(2,1)));

    figure; pie(fractions(:,1),{'',''})
    
    set(gcf,'Units','centimeters')
    set(gcf,'paperUnits','centimeters')

    pos = get(gcf,'Position');

    pos(3) = 17.8/2;
    pos(4) = 23.0/4;

    set(gcf,'Position',pos);
    set(gcf,'paperPosition',pos);
    pos_ax2 = [0.025, 0.015, 0.62,0.9];
    set(gca,'Position',pos_ax2)
    
    %For simulated Data:
    text(-1.2,1.15,[num2str(round(fractions(1,1)*100,1)), '\pm', num2str(round(fractions(1,2)*100,1)),'%'],'FontSize',10)


    text(1.1,-0.25,[num2str(round(fractions(2,1)*100,1)), '\pm', num2str(round(fractions(2,2)*100,1)),'%'],'FontSize',10)
    text(0.9,-0.5,['\tau = ',num2str(round(resTimes(1,1),3)), '\pm', num2str(round(resTimes(1,2),3)),' s'],'FontSize',10)
end
cmap(1,:) = [0.9,0.9,0.9];
if ~isempty(varargin)
    minRes = varargin{1}(1);
    maxRes = varargin{1}(2);
    
    cmap_base = colormap('jet');
    
    
    for i = 1:size(resTimes,1)
        ratio1 = (resTimes(i,1) - minRes)/(maxRes - minRes);
        ratio1 = max(ratio1,0);
        ratio1 = min(ratio1,1);
        indx = round(ratio1*63) + 1;
        cmap(i+1,:) = cmap_base(indx,:);
    end
    
else
    cmap(2,:) = [0,1,0];
    if size(resTimes,1) == 2
        cmap(3,:) = [1,1,0];
    end
end
colormap(cmap);