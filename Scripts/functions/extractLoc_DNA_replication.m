function [cellSummary]=extractLoc_DNA_replication(cellSummary,mskCell,propROI,imgH,imgW,pixelSize)
%% Get localisation of SPP1 DNA replication area from the nearest cell pole.
%% Results are written in cellSummary array (columns 12 to 15).

disp('DNA replication localization ...');
nROI=size(cellSummary,1);

for iROI=1:nROI
    if (cellSummary(iROI,3)==1)
        [ii,ij]=find(mskCell==iROI); % pixel of current cell 
        %img_profile=mskCell;img_profile(mskCell~=iROI)=NaN;
        
        % define profile line corresponding to the medial axis of bacteria
        % (and cut to the poles)
        angleROI=-propROI(iROI).Orientation; % cell orientation (axis image)
        iC=mean(ii);jC=mean(ij);
        lCell=sqrt((min(ii)-max(ii))^2+(min(ij)-max(ij))^2); % cell length
        %hold on; plot(jC,iC,'wo');
                
        i0=iC-0.6*lCell*sind(angleROI);% could we test is line is still in cur ROI to have a complete profile?
        j0=jC-0.6*lCell*cosd(angleROI);
        i1=iC+0.6*lCell*sind(angleROI);
        j1=jC+0.6*lCell*cosd(angleROI);
        figure(16);hold on; plot([j0,j1],[i0,i1],'w-.');
        jProfile=j0:sign(j1-j0):j1;
        iProfile=i0:sign(i1-i0):i1;
        
        if(numel(iProfile)>numel(jProfile))
            jProfile=linspace(jProfile(1),jProfile(end),numel(iProfile));
        else
            iProfile=linspace(iProfile(1),iProfile(end),numel(jProfile));
        end
        jProfile=floor(0.5+jProfile);
        iProfile=floor(0.5+iProfile);
        
        jProfile((iProfile>imgH) | (iProfile<1))=[];
        iProfile((iProfile>imgH) | (iProfile<1))=[];
        iProfile((jProfile>imgW) | (jProfile<1))=[];
        jProfile((jProfile>imgW) | (jProfile<1))=[];
        %hold on; plot(jProfile,iProfile,'w.');
        
        lProfile=sqrt((iProfile(1)-iProfile(end))^2 + (jProfile(1)-jProfile(end))^2)*pixelSize;
        nP=numel(iProfile);
        intCell=zeros(nP,3);
        intCell(:,1)=linspace(0,lProfile,nP);
        intCell(:,3)=1:nP;
        %intGP=intSPP1;
        for iP=1:nP
            intCell(iP,2)=mskCell(iProfile(iP),jProfile(iP));
            %   intGP(iP,2)=imgYFP_profile(iProfile(iP),jProfile(iP));
        end%for
        intCell(intCell(:,2)~=iROI,:)=[]; % remove part of the profile line outside of the cell
        
        % get nearest pole of the center of SPP1 DNA replication area
        ijProfile=[iProfile(intCell(:,3));jProfile(intCell(:,3))];        
        infCM=cellSummary(iROI,[7,6]);        
        poleDist(1,1)=sqrt((ijProfile(1,1)-infCM(1))^2+(ijProfile(2,1)-infCM(2))^2);
        poleDist(2,1)=sqrt((ijProfile(1,end)-infCM(1))^2+(ijProfile(2,end)-infCM(2))^2);
        
        % update cellSummary array with minimal distance pole-center, pole
        % positions and cell length
        cellSummary(iROI,12)=min(poleDist);
        if (find(poleDist==min(poleDist))==1)
            cellSummary(iROI,14:15)=[ijProfile(2,1),ijProfile(1,1)];
        else
            cellSummary(iROI,14:15)=[ijProfile(2,end),ijProfile(1,end)];
        end
        cellSummary(iROI,13)=(intCell(end,1)-intCell(1,1))/pixelSize;
    end
end
disp('DNA replication localization done');
end