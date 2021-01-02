
function JDH_PhB = correctBleachJDH(JDH, BleachRate, tlist,rlist, plotFlag)
% correctBleachJDH
% This function corrects a Jump distance histogram for photobleaching
% described by a single exponential.
% Input:
% The histogram JDH has the following structure
%
%           |  JDH(r1, t1)      JDH(r2, t1)    ...      JDH(rmax, t1)   |    
%           |  JDH(r1, t2)      JDH(r2, t2)    ...      JDH(rmax, t2)   |
%  JDH    = |      ...              ...        ...           ...        |
%           |  JDH(r1, tmax)    JDH(r2, tmax)   ...     JDH(rmax, tmax) |
%
%--------------------------------------------------------------------------
% tlist is a column vector with the times [t1, ...,tmax]
% Bleach rate is a scalar (units s^-1)


Rate1 = BleachRate(1);
Rate2 = BleachRate(2);
f1 = BleachRate(3);

CorrMatrix = repmat(tlist, [1 length(JDH(1,:))]);
CorrMatrix = f1*exp(-CorrMatrix*Rate1) + (1-f1)*exp(-CorrMatrix*Rate2);
JDH_PhB = JDH./CorrMatrix;

if plotFlag == 1
figure;
    imagesc(rlist, tlist, JDH_PhB);
    set(gca, 'FontSize', 12);
    
    xlabel('Jump distance [\mum]');
    ylabel('Time [s]');
    title('Jump Distance Histogram Distribution [Photobleached corrected]]');
    colorbar('location','EastOutside','XTick',0.5,'XTickLabel','Counts');
    
    axis image
end
end
