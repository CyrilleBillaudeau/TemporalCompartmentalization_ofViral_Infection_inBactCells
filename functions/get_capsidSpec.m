function [locCapsids]=get_capsidSpec(mskLocCapsids,mskCell,cellSummary,imgYFP,imgInfectionGlobal,imgH,imgW,nROI)
mskLocCapsids=bwlabel(mskLocCapsids>0);
propCapsids=regionprops(mskLocCapsids>0,'Area','Centroid','PixelList','MajorAxisLength','MinorAxisLength','Perimeter');
nLocCapsids=numel(propCapsids);

locCapsids=zeros(nLocCapsids,16);
shapeCapsids=zeros(nLocCapsids,5);

for iL=1:nLocCapsids
    %shapeCapsids(iL,:)=[bwarea(mskLocCapsids==iL),propCapsids(iL).Perimeter,-1,propCapsids(iL).MajorAxisLength,propCapsids(iL).MinorAxisLength];
    shapeCapsids(iL,:)=[propCapsids(iL).Area,propCapsids(iL).Perimeter,-1,propCapsids(iL).MajorAxisLength,propCapsids(iL).MinorAxisLength];
end%for
shapeCapsids(:,3)=shapeCapsids(:,5)./shapeCapsids(:,4);

for iL=1:nLocCapsids
    posCaps=round(propCapsids(iL).Centroid);
    msk_rows=posCaps(2)-3:posCaps(2)+3;msk_rows(msk_rows<1)=[];msk_rows(msk_rows>size(mskCell,1))=[];
    msk_cols=posCaps(1)-3:posCaps(1)+3;msk_cols(msk_cols<1)=[];msk_cols(msk_cols>size(mskCell,2))=[];
    curROI=mskCell(msk_rows, msk_cols);curROI=curROI(:);curROI(curROI==0)=[];
    curROI=median(curROI);
    if (mod(curROI,1)~=0)
        curROI=[];
    end
    
    if ((numel(curROI==1)) & (cellSummary(curROI,3)==1))
        
        locCapsids(iL,1:2)=[iL,curROI];
        
        % mass center of capsid and mean intensity (cols 3-5)
        [pixROI]=propCapsids(iL).PixelList;
        nPix=size(pixROI,1);
        for iPix=1:nPix; pixROI(iPix,3)=imgYFP(pixROI(iPix,2),pixROI(iPix,1));end
        xCM=sum(pixROI(:,1).*pixROI(:,3))/sum(pixROI(:,3));
        yCM=sum(pixROI(:,2).*pixROI(:,3))/sum(pixROI(:,3));
        locCapsids(iL,3:4)=[xCM,yCM];
        locCapsids(iL,5)=mean(pixROI(:,3));
        
        % distance dna replication - capsid assembly (col 6)
        [xyCMdna]=cellSummary(cellSummary(:,1)==locCapsids(iL,2),6:7);
        locCapsids(iL,6)=sqrt((xCM-xyCMdna(1))^2+(yCM-xyCMdna(2))^2);
        
        % overlap btw dna & capsids
        iCell=locCapsids(iL,2);
        
        curCell_DNA=(mskCell==iCell) & (imgInfectionGlobal ==3);%figure(205);clf;imagesc(curCell_DNA);axis equal;colormap(gray);
        curCell_Capsid=mskLocCapsids==iL;%figure(205);clf;imagesc(curCell_Capsid);axis equal;colormap(gray);
        curCell_Overlay=curCell_Capsid & curCell_DNA;%figure(205);clf;imagesc(curCell_Overlay);axis equal;colormap(gray);
        
        [iiDNA,~]=find(curCell_DNA>0);
        [iiCaps,~]=find(curCell_Capsid>0);
        [iiOverlay,~]=find(curCell_Overlay>0);
        
        sDNA=size(iiDNA,1);
        sCaps=size(iiCaps,1);
        sOverlay=size(iiOverlay,1);
        if (sOverlay==0)
            locCapsids(iL,7)=0;
            locCapsids(iL,8)=0;
        else
            locCapsids(iL,7)=100*sOverlay/(sCaps);
            locCapsids(iL,8)=100*sOverlay/(sDNA);
        end
        %            kColor=[0 0 0]/255;
        %            bColor=[0 130 254]/255;
        oColor=[255 128 1]/255;
        %            mColor=[255 0 255]/255;
        if (locCapsids(iL,7)>0)
            figure(201);hold on; text(locCapsids(iL,3)-0,locCapsids(iL,4)-10,num2str(locCapsids(iL,7),3),'Color',oColor,'FontName','Arial','FontWeight','bold');
            figure(201);hold on; plot(locCapsids(iL,3),locCapsids(iL,4),'k+');
        else
            figure(201);hold on; text(locCapsids(iL,3)+0,locCapsids(iL,4)+10,num2str(locCapsids(iL,7),3),'Color','red','FontName','Arial','FontWeight','bold');
            figure(201);hold on; plot(locCapsids(iL,3),locCapsids(iL,4),'ko');
        end
        
        %%
        j0=xCM;j1=xyCMdna(1);i0=yCM;i1=xyCMdna(2);
        imDNA=bwlabel(imgInfectionGlobal==3);
        curIDdna=imDNA(round(xyCMdna(2)),round(xyCMdna(1)));
        
        imProfile=zeros(size(mskCell));
        imProfile(imDNA==curIDdna)=1;
        imProfile(mskLocCapsids==iL)=2;
        
        %             figure(21);clf;imagesc(imProfile);axis equal;colormap(gray);
        %             figure(21);hold on; plot([j0,j1],[i0,i1],'y-.');
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
        lProfile=sqrt((iProfile(1)-iProfile(end))^2 + (jProfile(1)-jProfile(end))^2);
        nP=numel(iProfile);
        intCell=zeros(nP,3);
        intCell(:,1)=linspace(0,lProfile,nP);
        intCell(:,3)=1:nP;
        %intGP=intSPP1;
        for iP=1:nP
            intCell(iP,2)=imProfile(iProfile(iP),jProfile(iP));
            %   intGP(iP,2)=imgYFP_profile(iProfile(iP),jProfile(iP));
        end%for
        intCell(intCell(:,2)~=0,:)=[];
        
        if (~isempty(intCell))
            locCapsids(iL,9)=intCell(end,1)-intCell(1,1);
        else
            locCapsids(iL,9)=0;
        end
        %%
        
        %
        xCaps=locCapsids(iL,3);
        yCaps=locCapsids(iL,4);
        xRep=cellSummary(locCapsids(iL,2),6);
        yRep=cellSummary(locCapsids(iL,2),7);
        xPole=cellSummary(locCapsids(iL,2),14);
        yPole=cellSummary(locCapsids(iL,2),15);
        
        
        xPole=xPole-xRep;
        xCaps=xCaps-xRep;
        if (xPole<0)
            xPole=-xPole;
            xCaps=-xCaps;
        end
        
        yPole=yPole-yRep;
        yCaps=yCaps-yRep;
        if (yPole<0)
            yPole=-yPole;
            yCaps=-yCaps;
        end
        
        %figure(900);clf;hold on;plot(0,0,'go');plot(xPole,yPole,'ks');plot(xCaps,yCaps,'r+');axis equal
        
        if (abs(xCaps)>abs(yCaps))
            testPos=xCaps;
        else
            testPos=yCaps;
        end
        if testPos>0;
            locCapsids(iL,16)=1;
        else
            locCapsids(iL,16)=2;
        end
        %title(num2str(locCapsids(iL,16)));pause()
        
    else
        
        locCapsids(iL,2:end)=-1;
    end
end%for

for iROI=1:nROI
    distZ=locCapsids(locCapsids(:,2)==iROI,6);
    [z,indS]=sort(distZ);
    %[(1:numel(indS))',indS]
    nS=numel(indS);inv_indS=zeros(nS,1);
    for iS=1:nS;inv_indS(indS(iS))=iS;end
    locCapsids(locCapsids(:,2)==iROI,10)=inv_indS;
end

locCapsids(:,11:15)=shapeCapsids;
end%function