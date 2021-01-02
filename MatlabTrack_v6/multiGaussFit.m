function [gparam, gfit, pval_Ftest] = multiGaussFit(data)

%generate a histogram from the data
MN = min(data);
MX = max(data);

[n,x] = hist(data,MN:MX);

%fit to single gaussian
coeff_LB = [0, 0, 0];
coeff_UB = [Inf, Inf, Inf];

fitfun = @(COEF) (gfit1_ls(COEF, x)) - (n);
COEF0 = [mean(data), std(data),length(data)];

options = optimset('FunValCheck','off', 'MaxIter', 1000, 'MaxFunEvals', 1000);

[gparam(1).par, resNorm, residuals,exitflag,output,lambda,jacobian] = lsqnonlin(fitfun,COEF0,coeff_LB,coeff_UB,options);

ci = nlparci(gparam(1).par,residuals,'jacobian',jacobian);
gparam(1).sigma = (ci(:,2) - ci(:,1))/2;
gparam(1).sigma = gparam(1).sigma';

gfit{1} = gfit1_ls(gparam(1).par,x);

%fit to sum of 2 gaussians
coeff_LB = [0, 0, 0, 0, 0, 0];
coeff_UB = [Inf, Inf, Inf, Inf, Inf, Inf];

fitfun = @(COEF) (gfit2_ls(COEF, x)) - (n);
COEF0 = [mean(data), std(data),0.5,2*mean(data),std(data),length(data)];

options = optimset('FunValCheck','off', 'MaxIter', 1000, 'MaxFunEvals', 1000);

[gparam(2).par, resNorm, residuals,exitflag,output,lambda,jacobian] = lsqnonlin(fitfun,COEF0,coeff_LB,coeff_UB,options);

ci = nlparci(gparam(2).par,residuals,'jacobian',jacobian);
gparam(2).sigma = (ci(:,2) - ci(:,1))/2;
gparam(2).sigma = gparam(2).sigma';

gfit{2}(1,:) = gfit2_ls(gparam(2).par,x);
gfit{2}(2,:) = gfit1_ls([gparam(2).par(1:2),gparam(2).par(3)*gparam(2).par(6)],x);
gfit{2}(3,:) = gfit1_ls([gparam(2).par(4:5), (1 - gparam(2).par(3))*gparam(2).par(6)],x);


%see if this is an improvement over 1 gaussian, and if it is keep going
[h_Ftest, pval_Ftest(1), fstat] = FtestModelCompare(n,gfit{1}(1,:),gfit{2}(1,:),3,6);

if h_Ftest == 1
    %fit to sum of 2 gaussians
    coeff_LB = [0, 0, 0, 0, 0, 0, 0, 0, 0];
    coeff_UB = [Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf];
    
    fitfun = @(COEF) (gfit3_ls(COEF, x)) - (n);
    COEF0 = [mean(data), std(data),0.3,2*mean(data),std(data),0.3,3*mean(data),std(data),length(data)];
    
    options = optimset('FunValCheck','off', 'MaxIter', 1000, 'MaxFunEvals', 1000);
    
    [gparam(3).par, resNorm, residuals,exitflag,output,lambda,jacobian] = lsqnonlin(fitfun,COEF0,coeff_LB,coeff_UB,options);
    
    ci = nlparci(gparam(3).par,residuals,'jacobian',jacobian);
    gparam(3).sigma = (ci(:,2) - ci(:,1))/2;
    gparam(3).sigma = gparam(3).sigma';
    
    gfit{3}(1,:) = gfit3_ls(gparam(3).par,x);
    gfit{3}(2,:) = gfit1_ls([gparam(3).par(1:2),gparam(3).par(3)*gparam(3).par(9)],x);
    gfit{3}(3,:) = gfit1_ls([gparam(3).par(4:5), gparam(3).par(6)*gparam(3).par(9)],x);
   
    gfit{3}(4,:) = gfit1_ls([gparam(3).par(7:8), (1 - gparam(3).par(3) - gparam(3).par(6))*gparam(3).par(9)],x);
    
    %see if this is an improvement over 2 gaussian, and if it is keep going
    [h_Ftest, pval_Ftest(2), fstat] = FtestModelCompare(n,gfit{2}(1,:),gfit{3}(1,:),6,9);
    
    if h_Ftest == 1
        %fit to sum of 2 gaussians
        coeff_LB = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        coeff_UB = [Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf];
        
        fitfun = @(COEF) (gfit4_ls(COEF, x)) - (n);
        COEF0 = [mean(data), std(data),0.25,2*mean(data),std(data),0.25,3*mean(data),std(data),0.25,4*mean(data),std(data),length(data)];
        
        options = optimset('FunValCheck','off', 'MaxIter', 1000, 'MaxFunEvals', 1000);
        
        [gparam(4).par, resNorm, residuals,exitflag,output,lambda,jacobian] = lsqnonlin(fitfun,COEF0,coeff_LB,coeff_UB,options);
        
        ci = nlparci(gparam(4).par,residuals,'jacobian',jacobian);
        gparam(4).sigma = (ci(:,2) - ci(:,1))/2;
        gparam(4).sigma = gparam(4).sigma';
        
        gfit{4}(1,:) = gfit4_ls(gparam(4).par,x);
        gfit{4}(2,:) = gfit1_ls([gparam(4).par(1:2),gparam(4).par(3)*gparam(4).par(12)],x);
        gfit{4}(3,:) = gfit1_ls([gparam(4).par(4:5), gparam(4).par(6)*gparam(4).par(12)],x);
        gfit{4}(4,:) = gfit1_ls([gparam(4).par(7:8), gparam(4).par(9)*gparam(4).par(12)],x);
        gfit{4}(5,:) = gfit1_ls([gparam(4).par(10:11), (1 - gparam(4).par(3) - gparam(4).par(6) - gparam(4).par(9))*gparam(4).par(9)],x);
        
        %see if this is an improvement over 1 gaussian, and if it is keep going
        [h_Ftest, pval_Ftest(3), fstat] = FtestModelCompare(n,gfit{3}(1,:),gfit{4}(1,:),9,12);
    end
end



function yhat = gfit1_ls(coeff,xlist)

yhat = coeff(3)*normpdf(xlist,coeff(1),coeff(2));

function yhat = gfit2_ls(coeff,xlist)

yhat = coeff(6)*(coeff(3)*normpdf(xlist,coeff(1),coeff(2)) + (1-coeff(3))*normpdf(xlist,coeff(4),coeff(5)));

function yhat = gfit3_ls(coeff,xlist)

yhat = coeff(9)*(coeff(3)*normpdf(xlist,coeff(1),coeff(2)) + coeff(6)*normpdf(xlist,coeff(4),coeff(5)) + (1 - coeff(3) - coeff(6))*normpdf(xlist,coeff(7),coeff(8)));


function yhat = gfit4_ls(coeff,xlist)

yhat = coeff(12)*(coeff(3)*normpdf(xlist,coeff(1),coeff(2)) + coeff(6)*normpdf(xlist,coeff(4),coeff(5)) + coeff(9)*normpdf(xlist,coeff(7),coeff(8)) + ...
    (1 - coeff(3) - coeff(6) - coeff(9))*normpdf(xlist,coeff(10),coeff(11)));
