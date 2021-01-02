

%%%Histones and Photobleaching Analysis
close all
%%Import data and Background Photobleaching
period=0.2;  %Period of Data acquisition in seconds 
st=4;
%Loading of Histograms from the tracking file and previous photobleaching
%kinetics
A=ResTimeHist;  
%B=CumNParticles;

zz1=[];
zz2=[];
for i=1:length(A(:,2))
    zz1=[zz1;repmat(A(i,1),A(i,2),1)];
end



[CDF,t,CDFup,CDFlo]=ecdf(zz1,'function','survivor','alpha',0.01,'bounds','on');

init=min(find(t>=st*period));   %Initial value corresponds to the max(ShortesTrack,Nmin)   
init1=init;

%%%%Fit 1 two exponential            

%% Triple exponential function of the form (f1*(exp(-beta*x))+f2*(exp(-gamma*x)+f3*(exp(-eta*x)))) with parameters f1=p(1) f2=p(2) f3=p(3) beta=p(4) gamma=p(5) eta=p(6) 
ft = @(p,xdata)(p(1)*(exp(-p(4).*xdata))+p(2)*(exp(-p(5).*xdata)+p(3)*(exp(-p(6).*xdata))));
%% Double exponential function of the form ((f1*(exp(-gamma*x)+f2*(exp(-eta*x)))) with parameters f1=p(1) f2=p(2) gamma=p(3) eta=p(4)
ft2= @(p,xdata)(p(1)*(exp(-p(3).*xdata)+p(2)*(exp(-p(4).*xdata))));

%%Fit it to the difference between predictive values and actual values
fitft = @(p)ft(p,t(init1:end))-CDF(init1:end);

%%Fit it to the difference between predictive values and actual values
fitft2 = @(p)ft2(p,t(init1:end))-CDF(init1:end);



%%Now the fit for the triple exponential
lb=[0,0,0,0,0,0];
ub=[Inf,Inf,Inf,Inf,Inf,Inf];
p0 = 0.1*ones(1,6); 
problem1 = createOptimProblem('lsqnonlin','x0',p0,'objective',fitft,...
    'lb',lb,'ub',ub);
% ms = MultiStart('PlotFcns',@gsplotbestf);
[Parameter1,errormulti1] = run(MultiStart,problem1,50);

options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
[~,~,Res,~,~,~,J] = lsqnonlin(fitft,p0,lb,ub,options);
CI_Parameter1=nlparci(Parameter1,Res,'jacobian',J);


%%Now the fit for the double exponential
lb=[0,0,0,0];
ub=[Inf,Inf,Inf,Inf];
p0 = 0.1*ones(1,4); 
problem2 = createOptimProblem('lsqnonlin','x0',p0,'objective',fitft2,...
    'lb',lb,'ub',ub);
% ms = MultiStart('PlotFcns',@gsplotbestf);
[Parameter2,errormulti2] = run(MultiStart,problem2,50);


[~,~,Res,~,~,~,J] = lsqnonlin(fitft2,Parameter2,lb,ub);
CI_Parameter2=nlparci(Parameter2,Res,'jacobian',J);
% final=find(CDF<0.01,1);
figure;
hax=axes; 
scatter(t,CDF,'.')
hold on
plot(t,CDFup,'b--')
plot(t,CDFlo,'b--')
plot(t,ft(Parameter1,t),'g')
set(gca,'xscale','log')
set(gca,'yscale','log')
title('Photobleaching Kinetics - 3 exponentials')
xlabel('Time (Secs)')
ylabel('Survival Distribution')
hold off

figure;
hax=axes; 
scatter(t,CDF,'.')
hold on
plot(t,CDFup,'b--')
plot(t,CDFlo,'b--')
hold on
plot(t,ft2(Parameter2,t),'g')
set(gca,'xscale','log')
set(gca,'yscale','log')
title('Photobleaching Kinetics - 2 exponentials')
xlabel('Time (Secs)')
ylabel('Survival Distribution')
hold off

run('PhotobleachingModelSelection_Evidence.m')
