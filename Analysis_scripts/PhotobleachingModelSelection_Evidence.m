%Define the probability Function for the double and triple exponential
%respectively

t1=t(init+1); %remove first point. 
t2=t(end);

ft = @(p,xdata)(p(1)*(p(4)*exp(-p(4).*xdata))+p(2)*(p(5)*exp(-p(5).*xdata)+p(3)*(p(6)*exp(-p(6).*xdata))));

ft2= @(p,xdata)(p(1)*(p(3)*exp(-p(3).*xdata)+p(2)*(p(4)*exp(-p(4).*xdata))));


C=1;

H1=@(x) C*ft(Parameter1,x);

%Renormalization of the pdf based on the experiment
C=integral(H1,t1,t2);
H1=@(x) ft(Parameter1,x)/C;

C2=1;

H2=@(x) C2*ft2(Parameter2,x);

%Renormalization of the pdf based on the experiment
C2=integral(H2,t1,t2);
H2=@(x) ft2(Parameter2,x)/C2;

Data=zz1(zz1>=t1);
P1=0;
P2=0;

J=unique(Data);

BICL1=0;
BICL2=0;

for i=1:length(J)
    
    siz=length(find(Data==J(i)));
    P1=P1+siz*log10(integral(H1,J(i)-period/2,J(i)+period/2));
    P2=P2+siz*log10(integral(H2,J(i)-period/2,J(i)+period/2));
    
    BICL1=siz*log(integral(H1,J(i)-period/2,J(i)+period/2))+BICL1;
    BICL2=siz*log(integral(H2,J(i)-period/2,J(i)+period/2))+BICL2;
end

BIC1=5*log(length(Data))-2*BICL1;  %Triple Exponential
BIC2=3*log(length(Data))-2*BICL2; %Doubple Exponential





sigma=(period/2)/length(Data);

error1=((ft(Parameter1,t(init:end))-CDF(init:end)).^2);
error2=((ft2(Parameter2,t(init:end))-CDF(init:end)).^2);

RSS1=sum(((ft(Parameter1,t(init:end))-CDF(init:end)).^2));
RSS2=sum(((ft2(Parameter2,t(init:end))-CDF(init:end)).^2));

sigma1=var(error1);
sigma2=var(error2);


BIC11=length(t(init:end))*log(RSS1/length(t(init:end)))+3*log(length(t(init:end)));
BIC12=length(t(init:end))*log(RSS2/length(t(init:end)))+2*log(length(t(init:end)));

%BIC11=abs(BIC11);
%BIC12=abs(BIC12);

Evid=10*(P1-P2);


formatSpec = 'The evidence for a triple exponential fit is %f Db\n';
fprintf(formatSpec,Evid)

formatSpec = 'The BIC1 for a triple exponential and a double exponential fit respectively is %f, %f \n';
fprintf(formatSpec,[BIC1,BIC2])









