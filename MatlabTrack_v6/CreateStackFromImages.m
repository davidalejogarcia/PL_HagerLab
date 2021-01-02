Fold = uigetdir(pwd,'Select Folder containing files to combine into image stack');
curdir = pwd;

if Fold ~= 0
    cd(Fold);
    ims = dir('*.tif');
    if isempty(ims)
        fprintf('Specified folder contains no images');
    else
        defaults = {'473','473_Stack', '561','561_Stack'};
        prompt = {'Channel 1 String:','Channel 1 Output Name:','Channel 2 String:','Channel 2 Output Name:'};
        dlgtitle = 'Strings to match for the individual channels';
        FileParam = inputdlg(prompt,dlgtitle, 1, defaults);
        tmpfold = matlabpath;
        sep = strfind(tmpfold,pathsep);
        tmpfold = tmpfold(1:sep-1);
        cd(tmpfold)
        TmpExists = dir('ImagesTmp');
        if isempty(TmpExists)
            mkdir(tmpfold,'ImagesTmp');
        end
        tmpfold = [tmpfold, filesep,'ImagesTmp'];
        cd(Fold);
        im_names1 = {};
        im_names2 = {};
        for i = 1:length(ims)
            
            if strfind(ims(i).name,FileParam{1})
                im_names1 = [im_names1;ims(i).name];
            else
                im_names2 = [im_names2;ims(i).name];
            end
        end
    end
    imFile1_tmp = [tmpfold,filesep,FileParam{2},'.tif'];
    imFile1 = [Fold,filesep,FileParam{2},'.tif'];
    imFile2_tmp = [tmpfold,filesep,FileParam{4},'.tif'];
    imFile2 = [Fold,filesep,FileParam{4},'.tif'];
    
    for i = 1:size(im_names1,1)
        I = imread(fullfile(Fold,im_names1{i,:}),'tif');
        imwrite(I,imFile1_tmp,'WriteMode','append');
    end
    copyfile(imFile1_tmp,imFile1);
    
    for i = 1:size(im_names2,1)
        I = imread(fullfile(Fold,im_names2{i,:}),'tif');
        imwrite(I,imFile2_tmp,'WriteMode','append');
    end
    copyfile(imFile2_tmp,imFile2);
    delete(imFile1_tmp,imFile2_tmp);
    fprintf('Done!\n');            
  
end