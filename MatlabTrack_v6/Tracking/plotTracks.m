function plotTracks(Tracks, imageIx)
    delete(findobj(gca,'Color','g'));
    ImIx = find (Tracks(:,3) == imageIx); % find tracks with a segment in the image
    if ~isempty(ImIx)
        pIx = Tracks(ImIx,4);           %find the corrisponding particle index;
        for i =pIx'
        plotIx = find(Tracks(:,4) == i);% & Tracks(:,3) <= imageIx);
        hold on;
        plot(Tracks(plotIx,1),Tracks(plotIx,2),'g','LineWidth',2);
        hold off;
        end;
    end