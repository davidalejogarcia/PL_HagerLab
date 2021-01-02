function [pparam, pfit, pval_Ftest] = multiPoissFit(data)

%generate a histogram from the data
MN = min(data);
MX = max(data);

[n,x] = hist(data,MN:MX);

%fit to single gaussian
coeff_LB = [0, 0];
coeff_UB = [Inf, Inf];

fitfun = @(COEF) (pfit1_ls(COEF, x)) - (n);
COEF0 = [mean(data),length(data)];

options = optimset('FunValCheck','off', 'MaxIter', 1000, 'MaxFunEvals', 1000);

[pparam(1).par, resNorm, residuals,exitflag,output,lambda,jacobian] = lsqnonlin(fitfun,COEF0,coeff_LB,coeff_UB,options);

ci = nlparci(pparam(1).par,residuals,'jacobian',jacobian);
pparam(1).sigma = (ci(:,2) - ci(:,1))/2;
pparam(1).sigma = pparam(1).sigma';

pfit{1} = pfit1_ls(pparam(1).par,x);

%fit to sum of 2 gaussians
coeff_LB = [0, 0, 0, 0];
coeff_UB = [Inf, Inf, Inf, Inf];

fitfun = @(COEF) (pfit2_ls(COEF, x)) - (n);
COEF0 = [mean(data),0.5,2*mean(data),length(data)];

options = optimset('FunValCheck','off', 'MaxIter', 1000, 'MaxFunEvals', 1000);

[pparam(2).par, resNorm, residuals,exitflag,output,lambda,jacobian] = lsqnonlin(fitfun,COEF0,coeff_LB,coeff_UB,options);

ci = nlparci(pparam(2).par,residuals,'jacobian',jacobian);
pparam(2).sigma = (ci(:,2) - ci(:,1))/2;
pparam(2).sigma = pparam(2).sigma';

pfit{2} = pfit2_ls(pparam(2).par,x);

%see if this is an improvement over 1 gaussian, and if it is keep going
[h_Ftest, pval_Ftest(1), fstat] = FtestModelCompare(n,pfit{1},pfit{2},2,4);

if h_Ftest == 1
    %fit to sum of 2 gaussians
    coeff_LB = [0, 0, 0, 0, 0, 0];
    coeff_UB = [Inf, Inf, Inf, Inf, Inf, Inf];
    
    fitfun = @(COEF) (pfit3_ls(COEF, x)) - (n);
    COEF0 = [mean(data),0.3,2*mean(data),0.3,3*mean(data),length(data)];
    
    options = optimset('FunValCheck','off', 'MaxIter', 1000, 'MaxFunEvals', 1000);
    
    [pparam(3).par, resNorm, residuals,exitflag,output,lambda,jacobian] = lsqnonlin(fitfun,COEF0,coeff_LB,coeff_UB,options);
    
    ci = nlparci(pparam(3).par,residuals,'jacobian',jacobian);
    pparam(3).sigma = (ci(:,2) - ci(:,1))/2;
    pparam(3).sigma = pparam(3).sigma';
    
    pfit{3} = pfit3_ls(pparam(3).par,x);
    
    %see if this is an improvement over 2 gaussian, and if it is keep going
    [h_Ftest, pval_Ftest(2), fstat] = FtestModelCompare(n,pfit{2},pfit{3},4,6);
    
    if h_Ftest == 1
        %fit to sum of 2 gaussians
        coeff_LB = [0, 0, 0, 0, 0, 0, 0, 0];
        coeff_UB = [Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf];
        
        fitfun = @(COEF) (pfit4_ls(COEF, x)) - (n);
        COEF0 = [mean(data), std(data),0.25,2*mean(data),std(data),0.25,3*mean(data),std(data),0.25,4*mean(data),std(data),length(data)];
        
        options = optimset('FunValCheck','off', 'MaxIter', 1000, 'MaxFunEvals', 1000);
        
        [pparam(4).par, resNorm, residuals,exitflag,output,lambda,jacobian] = lsqnonlin(fitfun,COEF0,coeff_LB,coeff_UB,options);
        
        ci = nlparci(pparam(4).par,residuals,'jacobian',jacobian);
        pparam(4).sigma = (ci(:,2) - ci(:,1))/2;
        pparam(4).sigma = pparam(4).sigma';
        
        pfit{4} = pfit4_ls(pparam(4).par,x);
        
        %see if this is an improvement over 1 gaussian, and if it is keep going
        [h_Ftest, pval_Ftest(3), fstat] = FtestModelCompare(n,pfit{3},pfit{4},6,8);
    end
end



function yhat = pfit1_ls(coeff,xlist)

yhat = coeff(2)*poisspdf(xlist,coeff(1));

function yhat = pfit2_ls(coeff,xlist)

yhat = coeff(4)*(coeff(2)*poisspdf(xlist,coeff(1)) + (1-coeff(2))*poisspdf(xlist,coeff(3)));

function yhat = pfit3_ls(coeff,xlist)

yhat = coeff(6)*(coeff(2)*poisspdf(xlist,coeff(1)) + coeff(4)*poisspdf(xlist,coeff(3)) + (1 - coeff(2) - coeff(4))*poisspdf(xlist,coeff(5)));


function yhat = pfit4_ls(coeff,xlist)

yhat = coeff(8)*(coeff(2)*poisspdf(xlist,coeff(1)) + coeff(4)*poisspdf(xlist,coeff(3)) + coeff(6)*poisspdf(xlist,coeff(5)) + ...
    (1 - coeff(2) - coeff(4) - coeff(6))*poisspdf(xlist,coeff(7)));
