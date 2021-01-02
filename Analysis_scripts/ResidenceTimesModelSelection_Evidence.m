t1=t(init+1); %remove first point. 
t2=t(final);


H1=@(xdata) Parameter1(4)*(Parameter1(1)*(Parameter1(2).*exp(-Parameter1(2)*xdata))+(1-Parameter1(1))*(Parameter1(3).*exp(-Parameter1(3)*xdata))).*exp(-PB.*xdata)+...
    Parameter1(4)*(Parameter1(1).*(exp(-Parameter1(2)*xdata))+(1-Parameter1(1)).*(exp(-Parameter1(3)*xdata)))*PB.*exp(-PB.*xdata);

C=integral(H1,t1-period/2,t2+period/2);
H1=@(xdata) (Parameter1(4)*(Parameter1(1)*(Parameter1(2).*exp(-Parameter1(2)*xdata))+(1-Parameter1(1))*(Parameter1(3).*exp(-Parameter1(3)*xdata))).*exp(-PB.*xdata)+...
    Parameter1(4)*(Parameter1(1).*(exp(-Parameter1(2)*xdata))+(1-Parameter1(1)).*(exp(-Parameter1(3)*xdata)))*PB.*exp(-PB.*xdata))/C;



H2=@(x) (Parameter2(1).*x.^(-Parameter2(2)))*PB.*exp(-PB*x)+(Parameter2(1)*Parameter2(2).*x.^(-Parameter2(2)-1)).*exp(-PB*x);
C2=integral(H2,t1-period/2,t2+period/2);

H2=@(x) ((Parameter2(1).*x.^(-Parameter2(2)))*PB.*exp(-PB*x)+(Parameter2(1)*Parameter2(2).*x.^(-Parameter2(2)-1)).*exp(-PB*x))/C2;

%Triple Exponential
H3 = @(xdata)(Parameter3(1)*((Parameter3(4)+PB)*exp(-(Parameter3(4)+PB).*xdata))+Parameter3(2)*((Parameter3(5)+PB)*exp(-(Parameter3(5)+PB).*xdata))+Parameter3(3)*((Parameter3(6)+PB)*exp(-(Parameter3(6)+PB).*xdata)));
C3=integral(H3,t1-period/2,t2+period/2);
H3 = @(xdata)(Parameter3(1)*((Parameter3(4)+PB)*exp(-(Parameter3(4)+PB).*xdata))+Parameter3(2)*((Parameter3(5)+PB)*exp(-(Parameter3(5)+PB).*xdata))+Parameter3(3)*((Parameter3(6)+PB)*exp(-(Parameter3(6)+PB).*xdata)))/C3;

Data=zz1(zz1>=t1);
P1=0;
P2=0;
P3=0;


BICL1=0;
BICL2=0;
BICL3=0;

J=unique(Data);

rep=zeros(length(J),1);
for i=1:length(J)
    
    siz=length(find(Data==J(i)));
    P1=P1+siz*log10(integral(H1,J(i)-period/2,J(i)+period/2));
    P2=P2+siz*log10(integral(H2,J(i)-period/2,J(i)+period/2));
    P3=P3+siz*log10(integral(H3,J(i)-period/2,J(i)+period/2));
    
    BICL1=siz*log(integral(H1,J(i)-period/2,J(i)+period/2))+BICL1;
    BICL2=siz*log(integral(H2,J(i)-period/2,J(i)+period/2))+BICL2;
    BICL3=siz*log(integral(H3,J(i)-period/2,J(i)+period/2))+BICL3;
end

P12=log10(1+10^(P2-P1))+P1+log10(1/2);
P23=log10(1+10^(P3-P2))+P2+log10(1/2);
P13=log10(1+10^(P3-P1))+P1+log10(1/2);

BIC1=3*log(length(Data))-2*BICL1;  %Double Exponential
BIC2=2*log(length(Data))-2*BICL2; %Power Law
BIC3=5*log(length(Data))-2*BICL3; %Triple Exponential

%Evidence in Db
Evid1=10*log10(1/3)+10*(P1-P23); 
Evid2=10*log10(1/3)+10*(P2-P13);
Evid3=10*log10(1/3)+10*(P3-P12);



sigma=(period/2)/length(Data);

error1=((ft(Parameter1,t(init:final))-CDF(init:final)).^2);
error2=((ft2(Parameter2,t(init:final))-CDF(init:final)).^2);
error3=((ft3(Parameter3,t(init:final))-CDF(init:final)).^2);

RSS1=sum(((ft(Parameter1,t(init:final))-CDF(init:final)).^2));
RSS2=sum(((ft2(Parameter2,t(init:final))-CDF(init:final)).^2));
RSS3=sum(((ft3(Parameter3,t(init:final))-CDF(init:final)).^2));

sigma1=var(error1);
sigma2=var(error2);
sigma3=var(error3);


BIC11=length(t(init:final))*log(RSS1/length(t(init:final)))+3*log(length(t(init:final)));
BIC12=length(t(init:final))*log(RSS2/length(t(init:final)))+2*log(length(t(init:final)));
BIC13=length(t(init:final))*log(RSS3/length(t(init:final)))+5*log(length(t(init:final)));


formatSpec = 'The evidence for a double exponential, power-law and a triple exponential fit is %f, %f, %f respectively\n';
fprintf(formatSpec,[Evid1,Evid2,Evid3])

formatSpec = 'The BIC1 for a double exponential, power-law and a triple exponential fit is %f, %f, %f respectively\n';
fprintf(formatSpec,[BIC1,BIC2,BIC3])



