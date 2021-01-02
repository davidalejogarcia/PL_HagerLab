function [P, C] = OneCmpExpSimPDF_CDF(t,k,N)


resTimes = exprnd(1/k,N,1);


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