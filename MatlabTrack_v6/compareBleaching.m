function compareBleaching

[process_files,pnames] = uigetfile('*.mat','Select the Photobleaching analysis that you want to compare','','MultiSelect','on');

if ~iscell(process_files)
    process_files = {process_files};
end
if process_files{1} ~= 0
    tags = cell(length(process_files),1);
    for i = 1:length(process_files)
        current_tag = inputdlg(['Input a tag for the file: ',process_files{i}],'Name for legend');
        tags{i,:} = current_tag{1,:};
        load(fullfile(pnames,process_files{i}));
        Npart{i} = Nparticles;
        Npart_fit{i} = fit_NPart;
        int_bgsub{i} = int_bg_sub_avg;
        int_fit{i} = fit_bgall;
    end
    
    figure; plot(Npart{1}(:,2));
    hold on
    for i = 2:length(process_files)
        plot(Npart{i}(:,2));
    end
%     legend(tags);
    title('Particle Number - raw')
    pos = get(gcf,'Position');
    pos2 = pos;
    pos2(1) = 300;
    set(gcf,'Position',pos2);
    
    figure; plot(Npart_fit{1}(:,2));
    hold on
    for i = 2:length(process_files)
        plot(Npart_fit{i}(:,2));
    end
    legend(tags);
    title('Particle Number - fit')
    pos2(1) = 900;
    set(gcf,'Position',pos2);
    
    figure; plot(int_bgsub{1}(:,1)./int_bgsub{1}(1,1));
    hold on
    for i = 2:length(process_files)
        plot(int_bgsub{i}(:,1)./int_bgsub{i}(1,1));
    end
    legend(tags);
    title('Intensity - raw')
    pos2(1) = 300;
    pos2(2) = 100;
    set(gcf,'Position',pos2);
    
    figure; plot(int_fit{1}(:,2)./int_fit{1}(1,2));
    hold on
    for i = 2:length(process_files)
        plot(int_fit{i}(:,2)./int_fit{i}(1,2));
    end
    legend(tags);
    title('Intensity - fit')
    pos2(1) = 900;
    set(gcf,'Position',pos2);
end
    