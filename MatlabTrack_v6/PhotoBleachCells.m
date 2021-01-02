function PhotoBleachCells

[process_files,pnames] = uigetfile('*.mat','Select the MatTrack files to use for photobleaching analysis','','MultiSelect','on');
% process_files1 = dir('*ROIs.mat');
% process_files2 = dir('*preprocess.mat');
% process_files = [process_files1;process_files2];
% process_files = process_files(3:2:length(process_files));
[outfile, outdir] = uiputfile('*.mat','Specify the name of the output file',fullfile(pnames,'PB_characteristics.mat'));

cell_int = cell(length(process_files),1);
bg_int = cell(length(process_files),1);
npx_cells_all = zeros(length(process_files),1);
npx_bg = zeros(length(process_files),1);

frameInterval = 0.2;


nImages = 0;
if ~iscell(process_files)
    if process_files == 0
        process_files = [];
    else
        process_files = {process_files};
    end
end
    

maxTpoints_Npart = 0;
if outfile ~= 0
    for i = 1:length(process_files)
        load(fullfile(pnames,process_files{i}));
        tVec = (0:Results.Data.nImages-1)';
        tVec = tVec*frameInterval;
        mask = Results.Process.ROIimage;
        mask_all = zeros(size(mask{1}));
        for j = 1:length(mask)
            mask_all(mask{j} > 0) = 1;
        end
        cell_int{i} = zeros(Results.Data.nImages,length(mask));
        cell_int_all{i} = zeros(Results.Data.nImages,1);
        bg_int{i} = zeros(Results.Data.nImages,1);


        for j = 1:Results.Data.nImages
            cell_int_all{i}(j,:) = sum(double(Results.Data.imageStack(j).data(mask_all == 1)));

            bg_int{i}(j,:) = sum(double(Results.Data.imageStack(j).data(mask_all == 0)));
            for k = 1:length(mask)
                cell_int{i}(j,k) = sum(double(Results.Data.imageStack(j).data(mask{k} > 0)));
            end
        end
        for j = 1:length(mask)
            npx_cells(i,j) = numel(mask{j}(mask{j} > 0));
        end
        npx_cells_all(i,:) = numel(mask_all(mask_all == 1));
        npx_bg(i,:) = numel(mask_all(mask_all == 0));
        nImages(i,:) =Results.Data.nImages;
        Nparticles{i}(:,1) = Results.PreAnalysis.NParticles(:,1);
        Nparticles{i}(:,2) = sum(Results.PreAnalysis.NParticles(:,2:end),2);
        
    end
    
    nParticles(:,1) = Nparticles{i}(:,1);
    nParticles(:,2) = zeros(size(nParticles(:,1),1),1);
    for i = 1:length(nParticles)
        for j = 1:length(Nparticles)
            nParticles(i,2) = nParticles(i,2) + Nparticles{j}(i,2);
        end
    end
    Nparticles = nParticles;
    int_bg_sub = zeros(max(nImages),length(process_files));
    for i = 1:length(process_files)
        int_bg_sub(1:nImages(i,:),i) = (cell_int_all{i}./npx_cells_all(i)) - (bg_int{i}./npx_bg(i));
    end
    int_bg_sub_avg = mean(int_bg_sub,2);

    int_bgmin_sub = zeros(max(nImages),length(process_files));
    for i = 1:length(process_files)
        int_bgmin_sub(1:nImages(i,:),i) = (cell_int_all{i}./npx_cells_all(i)) - min((bg_int{i}./npx_bg(i)));
    end
    int_bgmin_sub_avg = mean(int_bgmin_sub,2);

    [fitpar_bgall,espSigma_bgall, fit_bgall]= ExpDecay_2Cmp_fit([tVec(1:nImages(i,:)),int_bg_sub_avg(1:nImages(i,:))], [1 0.1]);
    [fitpar_bgmin,espSigma_bgmin, fit_bgmin]= ExpDecay_2Cmp_fit([tVec(1:nImages(i,:)),int_bgmin_sub_avg(1:nImages(i,:))], [1 0.1]);
    [fitpar_NPart,espSigma_NPart, fit_NPart]= ExpDecay_2Cmp_fit([Nparticles(:,1),Nparticles(:,2)], [1 0.1]);
    % save('PB_characteristics.mat','cell_int','bg_int','npx*','int_bg_sub*','int_bgmin_sub*','fit*','esp*','Nparticles');
    save(fullfile(outdir,outfile),'cell_int','cell_int_all','bg_int','npx*','int_bg_sub*','int_bgmin_sub*','fit*','esp*','Nparticles');
%     figure; plot(int_bg_sub_avg./int_bg_sub_avg(1))
%     hold on
%     plot(fit_bgall(:,2)./fit_bgall(1,2))
end