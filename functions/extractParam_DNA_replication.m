function cellSummary=extractParam_DNA_replication(cellSummary,propROI,imgRFP,imgInfectionGlobal)
%% Extract localization and morphological properties of SPP1 DNA replication
%% area for all mono-infected cells. 
%% Results are written in cellSummary array (columns 6 to 11).

disp('DNA replication area parameters measurments ...');
nROI=size(cellSummary,1);
for iROI=1:nROI
    if (cellSummary(iROI,3)==1) % mono-infected cells only        
        % Get SPP1 DNA replication area inside the current cell
        pixROI=propROI(iROI).PixelList;% cell area
        nPix=size(pixROI,1);
        pixROI=[pixROI,zeros(nPix,1)];        
        for iPix=1:nPix; pixROI(iPix,3)=imgInfectionGlobal(pixROI(iPix,2),pixROI(iPix,1));end
        pixROI(pixROI(:,3)<3,:)=[];
        nPix=size(pixROI,1);
        areaDNAReplicationFocus=nPix;% area (pixel unit)
        
        % Get the centroid of SPP1 DNA replication area (weighted by fluo.
        % intensity) (pixel unit)
        for iPix=1:nPix; pixROI(iPix,3)=imgRFP(pixROI(iPix,2),pixROI(iPix,1));end
        xCM=sum(pixROI(:,1).*pixROI(:,3))/sum(pixROI(:,3));
        yCM=sum(pixROI(:,2).*pixROI(:,3))/sum(pixROI(:,3));
        %hold on; plot(xCM,yCM,'m+');
        %text(posROI(1)-20,posROI(2)-20,num2str(iROI),'Color','w')
        
        % Get region properties of SPP1 DNA replication area (major and
        % minor axis, assuming an ellipse shape (pixel unit)
        pixROI(:,1)=pixROI(:,1)-min(pixROI(:,1))+1;
        pixROI(:,2)=pixROI(:,2)-min(pixROI(:,2))+1;
        curDNAReplication=zeros(max(pixROI(:,1)),max(pixROI(:,2)));
        for iPix=1:nPix; curDNAReplication(pixROI(iPix,1),pixROI(iPix,2))=1;end
        res=regionprops(curDNAReplication,'MajorAxisLength','MinorAxisLength');
        majorAxisDNAReplication=res.MajorAxisLength;
        minorAxisDNAReplication=res.MinorAxisLength;
        
        % update cellSummary array
        cellSummary(iROI,6:8)=[xCM,yCM,areaDNAReplicationFocus];     
        cellSummary(iROI,9:11)=[majorAxisDNAReplication,minorAxisDNAReplication,mean(pixROI(:,3))];
    end%if
end
disp('DNA replication area parameters measurments done');
end