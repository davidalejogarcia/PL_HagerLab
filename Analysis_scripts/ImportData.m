
clear
[file,path] = uigetfile('*.xlsx');
if isequal(file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(path,file)]);

   T = readtable(fullfile(path,file));
   ResTimeHist=T{:,:};
end

[~,name,~] = fileparts(fullfile(path,file));

save(name,'ResTimeHist')


clear all
clc