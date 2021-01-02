close all



%%Import data and Background Photobleaching
period=0.2;  %Period of Data acquisition in seconds 
st=5;
PB=0.078; 
%If the whole data is to be analyzed set final to 1. Otherwise, specify the
%minimum number of accumulated tracks
final=1;
%Loading of Histograms from the tracking file and previous photobleaching
%kinetics
A=ResTimeHist;  

Results.totalTracks=sum(A(:,2));
zz1=[];
zz2=[];
for i=1:length(A(:,2))
    zz1=[zz1;repmat(A(i,1),A(i,2),1)];
end
[CDF,t,CDFup,CDFlo]=ecdf(zz1,'function','survivor','alpha',0.01,'bounds','on');
st=find(t==st*period);
edge=period-period/2:period:max(zz1)+period;
h=histogram(zz1,'Normalization','cumcount','BinEdges',edge);
NormCount=h.Values(end)-h.Values;
if final==1
    final=length(t);
else
    final2=min(find(NormCount<=final));
    final=find(t==final2*period);
end
close all
%%Photobleaching Correction
P2=exp(-PB*t); %H2B Photobleaching parameter

CDF=CDF./P2;
CDFup=CDFup./P2;
CDFlo=CDFlo./P2;
st=st(end);

CDFlo=CDFlo/CDFlo(st);
CDFup=CDFup/CDFup(st);
CDF=CDF/CDF(st);


%Previous fit 
init=st;   %Initial value corresponds to the min(ShortesTrack,Nmin)   
            
%%%%Fit 1 two exponential            
%% Double exponential function of the form (A*(f*(exp(-gamma*x)+(1-f)*(exp(-eta*x)))) with parameters f=p(1) gamma=p(2) eta=p(3) A=p(4)
ft= @(p,xdata)p(4)*(p(1)*(exp(-p(2)*xdata))+(1-p(1))*(exp(-p(3)*xdata)));

%%Fit it to the difference between predictive values and actual values
fitft = @(p)ft(p,t(init:final))-CDF(init:final);

%%Now the fit for the double exponential
lb=[0,0,0,0];
ub=[1,Inf,Inf,Inf];
p0 = 0.1*ones(1,4); 

problem1 = createOptimProblem('lsqnonlin','x0',p0,'objective',fitft,...
    'lb',lb,'ub',ub);
[Parameter1,errormulti1] = run(MultiStart,problem1,50);

options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
[~,~,Res,~,~,~,J] = lsqnonlin(fitft,p0,lb,ub,options);
CI_Parameter1=nlparci(Parameter1,Res,'jacobian',J);


figure;
scatter(t(1:final),CDF(1:final),'.')
hold on
plot(t(1:final),CDFup(1:final),'b--')
plot(t(1:final),CDFlo(1:final),'b--')
hold on
plot(t,ft(Parameter1,t),'g')
set(gca,'xscale','log')
set(gca,'yscale','log')
title('Fit 1 (Double Exponential)')
xlabel('Time (Secs)')
ylabel('Survival Distribution')
hold off


%%%%Fit 2 Power-Law
        

%% Power Law Fit of the form A*x^(-b) with parameter A=p(1) b=p(2)
ft2= @(p,xdata)(p(1).*xdata.^(-p(2)));

%%Fit it to the difference between predictive values and actual values
fitft2 = @(p)ft2(p,t(init:final))-CDF(init:final);

%%Now the fit for the double exponential
lb=[0,0];
ub=[Inf,Inf];
p0 = 0.1*ones(1,2); 
problem2 = createOptimProblem('lsqnonlin','x0',p0,'objective',fitft2,...
    'lb',lb,'ub',ub);
[Parameter2,errormulti2] = run(MultiStart,problem2,50);

[~,~,Res,~,~,~,J] = lsqnonlin(fitft2,Parameter2,lb,ub);
CI_Parameter2=nlparci(Parameter2,Res,'jacobian',J);

figure;
scatter(t(1:final),CDF(1:final),'.')
hold on
plot(t(1:final),CDFup(1:final),'b--')
plot(t(1:final),CDFlo(1:final),'b--')
plot(t,ft2(Parameter2,t),'g')
set(gca,'xscale','log')
set(gca,'yscale','log')
title('Fit 2 (Power-Law)')
xlabel('Time (Secs)')
ylabel('CDF')
hold off


%%%Fit Triple Exponential

%% Triple exponential function of the form (f1*(exp(-beta*x))+f2*(exp(-gamma*x)+f3*(exp(-eta*x)))) with parameters f1=p(1) f2=p(2) f3=p(3) beta=p(4) gamma=p(5) eta=p(6) 
ft3 = @(p,xdata)(p(7)*(p(1)*(exp(-p(4).*xdata))+p(2)*(exp(-p(5).*xdata)+p(3)*(exp(-p(6).*xdata)))));
fitft3 = @(p)ft3(p,t(init:final))-CDF(init:final);


lb=[0,0,0,0,0,0,0];
ub=[Inf,Inf,Inf,Inf,Inf,Inf,Inf];
p0 = 0.1*ones(1,7); 
problem3 = createOptimProblem('lsqnonlin','x0',p0,'objective',fitft3,...
    'lb',lb,'ub',ub);
[Parameter3,errormulti3] = run(MultiStart,problem3,50);

options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
[~,~,Res,~,~,~,J] = lsqnonlin(fitft3,p0,lb,ub,options);
CI_Parameter3=nlparci(Parameter3,Res,'jacobian',J);

figure;
hax=axes; 
scatter(t(1:final),CDF(1:final),'.')
hold on
plot(t(1:final),CDFup(1:final),'b--')
plot(t(1:final),CDFlo(1:final),'b--')
plot(t,ft3(Parameter3,t),'g')
set(gca,'xscale','log')
set(gca,'yscale','log')
title('fit3 (3 exponentials)')
xlabel('Time (Secs)')
ylabel('Survival Distribution')
hold off


run('ResidenceTimesModelSelection_Evidence.m')

