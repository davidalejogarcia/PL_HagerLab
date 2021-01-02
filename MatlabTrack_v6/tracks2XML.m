function tracks2XML(Tracks,filename,snrval,denStr)

%Writes an XML document from the tracks identified with integratedTrackGui
%for use with the Particle Tracking Challege analysis software.

docNode = com.mathworks.xml.XMLUtils.createDocument('root');
toc = docNode.getDocumentElement;

product = docNode.createElement('TrackContestISBI2012');

product.setAttribute('snr',num2str(snrval));
product.setAttribute('density',denStr);

for i = 1:max(Tracks(:,4))
    Track_tmp = Tracks(Tracks(:,4) == i,:);
    if ~isempty(Track_tmp)
        part = docNode.createElement('particle');
        for j = 1:size(Track_tmp,1)
            det = docNode.createElement('detection');
            det.setAttribute('t',num2str(Track_tmp(j,3)-1));
            det.setAttribute('x',num2str(Track_tmp(j,1)-1));
            det.setAttribute('y',num2str(Track_tmp(j,2)-1))
            det.setAttribute('z','0');
            part.appendChild(det);
        end
        
        product.appendChild(part);
    end
end
toc.appendChild(product);
xmlwrite(filename,docNode);

