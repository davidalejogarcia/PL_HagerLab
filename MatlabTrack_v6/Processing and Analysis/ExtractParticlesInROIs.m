function Data = ExtractParticlesInROIs()

[fname, pname] = uigetfile('*.mat','Select the Mat Track files from which data should be extracted','MultiSelect','on');
[fname_out,pname_out] = uiputfile('*.txt','Save Extracted Particle Numbers','ParticlesInROI.txt');
if ~iscell(fname)
    fname = {fname};
end
Data.FileName = {};
Data.Nparticles = [];
if fname{1} ~= 0
    for i = 1:size(fname,2)
        clear nPart_tmp;
        Tmp = load(fullfile(pname,fname{i}));
        img_name = Tmp.Results.Data.fileName;
        if isfield(Tmp.Results.Tracking,'Particles')
            Particles = Tmp.Results.Tracking.Particles;
            Particles(Particles(:,1) == 0,:) = [];
            Particles = Particles(Particles(:,6) == 1,13);
        else
            Particles = Tmp.Results.Tracking.Centroids;
            Particles(Particles(:,1) == 0,:) = [];
            Particles = Particles(Particles(:,6) == 1,7);
        end
        for j = 1:max(Particles)
            nPart_tmp(j,1) = j;
            nPart_tmp(j,2) = size(Particles(Particles == j),1);
        end
        fileName_tmp = repmat({img_name},size(nPart_tmp,1),1);
        
        Data.FileName = [Data.FileName; fileName_tmp];
        Data.Nparticles = [Data.Nparticles; nPart_tmp];
        
        
        
    end
    
    
    if fname_out ~= 0
        fid = fopen(fullfile(pname_out,fname_out),'w');
        fprintf(fid,'Image Name\tROI Number\tParticles\n');
        for i = 1:size(Data.FileName,1)
            fprintf(fid,'%s\t%u\t%u\n',Data.FileName{i,:},Data.Nparticles(i,1),Data.Nparticles(i,2));
        end
        fclose(fid);
    end
    
end