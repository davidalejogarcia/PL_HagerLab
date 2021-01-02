%%Merging of Experiments
clear
clc
d=dir('*.mat');
eta=[];
 zz1=[];
for ix=1:length(d)
   
  fn=d(ix).name;
  if strcmp(fn(1),'.')==0
      load(fn)
      A=ResTimeHist;  
      B=CumNParticles;

        for i=1:length(A(:,2))
            zz1=[zz1;repmat(A(i,1),A(i,2),1)];
        end 
  end
      
end

[CDF,t,CDFup,CDFlo]=ecdf(zz1,'function','survivor','alpha',0.01,'bounds','on');

Temp(:,1)=t;
Temp(:,2)=CDF*length(zz1);
CumNParticles=Temp; %Approximation

[a,b]=hist(zz1,unique(zz1));
Temp2(:,1)=b;
Temp2(:,2)=a.';
ResTimeHist=Temp2;


save('Merged.mat','CumNParticles','ResTimeHist')
close all 
clear 
clc
