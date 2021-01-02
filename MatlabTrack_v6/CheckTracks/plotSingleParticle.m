function plotSingleParticle (Track, imageIx, axesHandle)

axes(axesHandle);
delete(findobj(gca,'Color','r'));

idx = find (Track(:,3) == imageIx);

if ~isempty(idx)
    hold on;
    plot(Track(idx,1), Track(idx,2),'or', 'MarkerSize', 12);
    hold off;
end