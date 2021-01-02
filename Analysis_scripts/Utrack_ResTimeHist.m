clear
clc
period=0.2;
d=dir('*.mat');
j=1;
for ix=1:length(d)
   
  fn=d(ix).name;
  if strcmp(fn(1),'.')==0
    load(fn)
    B=[];
    B=struct2cell(tracksFinal);
    for i=1:length(B)
        Surv(j)=period*length(B{1,i});
        j=j+1;
    end
  end
      
end

A=unique(Surv);
N=histc(Surv,A);

ResTimeHist(:,1)=A.';
ResTimeHist(:,2)=N.';

save('UtrackMergedData','ResTimeHist')