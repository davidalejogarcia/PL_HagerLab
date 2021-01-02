
%open dialog to select folders
[fnames, pname] = uigetfile('*.mat','Open MatTrack files containing ROI to extract','Open ROIs','MultiSelect','on');

if ~isempty(fnames)
    %load in the first filename to determine which ROI to use
    NoiseAddAns = questdlg('Do you want to add white noise to the resulting image?','Add noise','Yes','No','No');
    ROIChoice = [];
    nFiles = length(fnames);
     ROIidx = 0;
    %determine the size of the ROI in each file
    for i = 1:nFiles
       
        Temp = load([pname,fnames{i}],'Results');
        if isfield(Temp.Results.Process,'ROIlabel')
            if size(Temp.Results.Process.ROIlabel,1) > 1
                if ~isempty(ROIChoice) %if a named ROI was selected before check to see if an ROI with that name exists in this dataset
                    if ~strcmp(ROIChoice{1},'All')
                        for j = 1:length(ROIChoice)
                            for k = 1:length(Temp.Results.Process.ROIlabel)
                                if strcmpi(ROIChoice{j},Temp.Results.Process.ROIlabel{k})
                                    ROIidx = k;
                                    break
                                end
                            end
                        end
                    else
                        ROIidx = size(Temp.Results.Process.ROIlabel,1) + 1; % If all was selected before, just set the index to something greater than the number of ROIs
                    end
                end
                
                if ROIidx == 0 %if nothing has been found automatically, provide the user with a dialog box to select the ROI
                    
                    ROIstring = Temp.Results.Process.ROIlabel;
                    
                    ROIidx = ROIchooseDlg(ROIstring);
                    
                end
                
                
                
            end
            for j = 1:Temp.Results.Data.nImages
                PxMean(i,j) = mean(Temp.Results.Data.imageStack(j).data(:));
                PxStd(i,j) = std(double(Temp.Results.Data.imageStack(j).data(:)));
            end
        end
        curROI = Temp.Results.Process.ROIimage{ROIidx,:};
        [x_cur,y_cur] = find(curROI == 1);
        extent(i,:) = [(max(x_cur) - min(x_cur) + 7), (max(y_cur) - min(y_cur) + 7)];
        
        
        cent(i,:) = [round(mean(double(x_cur))),round(mean(double(y_cur)))];
%         store the image
        Stacks{i} = Temp.Results.Data.imageStack;
        nImages(i) = Temp.Results.Data.nImages;
        ROIimage{i} = curROI;
        
    end
    allMeans = mean(PxMean);
    allStds = mean(PxStd);
    savefold = [pname,'mergedROI_',Temp.Results.Process.ROIlabel{ROIidx,:}];
    if strcmp(NoiseAddAns,'Yes')
        savefold = [savefold, '_Noise'];
    end
        mkdir(savefold);
    ROIborders = round(max(extent)/2);
    maxT = min(nImages);
    squareSide = ceil(sqrt(nFiles));
    nPanels = squareSide.^2;
    for i = 1:maxT
        if strcmp(NoiseAddAns,'Yes')
            allROI = normrnd(allMeans(i),allStds(i),2*ROIborders(1)+1,2*ROIborders(2)+1,nPanels);
        else
            allROI = zeros(2*ROIborders(1)+1,2*ROIborders(2)+1,nPanels)+allMeans(i);
        end
        allROI2 = reshape(allROI,squareSide*(2*ROIborders(1)+1),squareSide*(2*ROIborders(2)+1));
%         allROI = [];
%         fillerROI = zeros(2*ROIborders(1)+1,2*ROIborders(2)+1);
        for j = 1:nFiles
            curROI = uint16(double(Stacks{j}(i).data).*double(ROIimage{j}));
            curROI_cut = curROI(cent(j,1)-ROIborders(1):cent(j,1)+ROIborders(1),cent(j,2)-ROIborders(2):cent(j,2)+ROIborders(2));
            zeroInd = find(curROI_cut == 0);
            if strcmp(NoiseAddAns,'Yes')
                for k = 1:length(zeroInd);
                    curROI_cut(zeroInd(k)) = normrnd(allMeans(i),allStds(i));
                end
            else
                curROI_cut(curROI_cut == 0) = allMeans(i);
            end
%             allROI(:,:,j) = curROI(cent(j,1)-ROIborders(1):cent(j,1)+ROIborders(1),cent(j,2)-ROIborders(2):cent(j,2)+ROIborders(2));
            allROI(:,:,j) = curROI_cut;
        end
%         for j = 1:nPanels-nFiles
%             allROI = [allROI; fillerROI];
%         end
%         allROI = reshape(allROI,squareSide*(2*ROIborders(1)+1),squareSide*(2*ROIborders(2)+1));
        ind = 1;
        for m = 1:squareSide
            for n = 1:squareSide
                allROI2((m-1)*size(allROI,1)+1:m*size(allROI,1),(n-1)*size(allROI,2)+1:n*size(allROI,2)) = allROI(:,:,ind);
                ind = ind+1;
            end
        end

        
        tpoint = num2str(i);
        if i < 10
            tpoint = ['0' tpoint];
        end
        if i < 100
            tpoint = ['0' tpoint];
        end
        if maxT >= 1000 && i < 1000
            tpoint = ['0' tpoint];
        end
        imwrite(uint16(allROI2),[savefold, filesep, 'mergedROI_t',tpoint,'.tif']);
    end
end