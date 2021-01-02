function [P, C,C2, N1,tvec, resTimes] = TwoCmpExpSimPDF_CDF(t,k1,k2,f1,N)

N1 = round(N*f1);
N2 = N - N1;

resTimes1 = exprnd(1/k1,N1,1);
resTimes2 = exprnd(1/k2,N2,1);
resTimes3 = [20;25; 30];

resTimes = [resTimes1;resTimes2];
dT = mean(diff(t));
tMin = min(min(t),min(resTimes));
tMax = max(max(t),max(resTimes));
P = hist(resTimes,tMin:dT:tMax);
tvec = tMin:dT:tMax;
if tMin < min(t)
    P(tvec < min(t)) = [];
end
if tMax > max(t)
    dEps = tvec - max(t);
    ind = find(dEps > 1000*eps,1,'first');
    P1 = P(1:ind-1);
    P1(end) = P1(end) + sum(P(ind:end));
    P = P1;
end

P = P./sum(P);
P = P/mean(diff(t));
P = P(:);

C = cumsum(P,'reverse');
C = C./C(1);
C = C(:);
C2 = [];
% 
% resTimes = [resTimes1;resTimes2;resTimes3];
% dT = mean(diff(t));
% tMin = min(min(t),min(resTimes));
% tMax = max(max(t),max(resTimes));
% P2 = hist(resTimes,tMin:dT:tMax);
% tvec2 = tMin:dT:tMax;
% if tMin < min(t)
%     P2(tvec2 < min(t)) = [];
% end
% if tMax > max(t)
%     dEps = tvec2 - max(t);
%     ind = find(dEps > 1000*eps,1,'first');
%     P12 = P2(1:ind-1);
%     P12(end) = P12(end) + sum(P2(ind:end));
%     P2 = P12;
% end
% 
% P2 = P2./sum(P2);
% P2 = P2/mean(diff(t));
% P2 = P2(:);
% 
% C2 = cumsum(P2,'reverse');
% C2 = C2./C2(1);
% C2 = C2(:);