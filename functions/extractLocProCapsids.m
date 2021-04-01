function [locCapsids]=extractLocProCapsids(cellSummary,mskCell,imgYFP,imgInfectionGlobal,minCapsArea,imgW,imgH,nK_caps)
%% Detect and localize procapsids in mono-infected cells only.
%%
%% Return a table with localization, size and morphological information of 
%% detected procapsids organized as follow:
%% col 1 - capsidID
%% col 2 - cellID
%% col 3 - capsid_xCM (weighted center)
%% col 4 - capsid_yCM (weighted center)
%% col 5 - capsid_avgInt
%% col 6 - capsid_replication_dist_CM
%% col 7 - capsid_replication_overlap1
%% col 8 - replication_capsid_overlap2
%% col 9 - capsid_replication_interdistance
%% col 10 - capsid_rank_dist_replication
%% col 11 - capsid_area
%% col 12 - capsid_perimeter
%% col 13 - capsid_aspectRatio
%% col 14 - capsid_majorAxisLength
%% col 15 - capsid_minorAxisLength
%% col 16 - zoneCapside_replicationPole

disp('Pro-Capsid loc parameters ...');

locCapsids=[];

if (~isempty(cellSummary))
    
    imgCapsids=imgYFP;%medfilt2(imgYFP,[3 3]);
    
    %% select only mono-infected cells    
    mskCellNotInfected=zeros(size(mskCell));
    nROI=size(cellSummary,1);
    for iROI=1:nROI
        if (cellSummary(iROI,3)==0)
            %disp(cellSummary(iROI,[1,3]))
            mskCellNotInfected(mskCell==iROI)=1;
        end
    end%for
    mskCellNotInfected=bwlabel(mskCellNotInfected);
    nROI_neg=max(mskCellNotInfected(:));
    intROI_neg=[];
    for iROI=1:nROI_neg
        int_neg=imgCapsids(mskCellNotInfected==iROI);
        intROI_neg=[intROI_neg;int_neg];
    end
    intROI_neg=mean(intROI_neg)+3*std(intROI_neg);
    
    mskCellMonoInfected=zeros(size(mskCell));
    nROI=size(cellSummary,1);
    for iROI=1:nROI
        if (cellSummary(iROI,3)==1)
            %disp(cellSummary(iROI,[1,3]))
            mskCellMonoInfected(mskCell==iROI)=1;    
        end
    end%for       
    %figure(20);clf;imagesc(mskCellMonoInfected);axis equal;colormap(gray);title('Mono Infected Cells')
    
    %% For viral procapsids, segmentation of procapsids foci was achieved by 
    % the k-means method (with 2 kernels), and by removing areas smaller 
    % than minCapsArea pixels² and areas displaying less than 5% pixels with intensity 
    % above an activation threshold. The activation threshold was defined 
    % on the fluorescence intensity of noninfected cells as the 
    % average + 3*standard deviation.
    
    % Gaussian filtering
    imgCapsids_Gauss=imgaussfilt(imgCapsids,2);
    figure(210);clf;imagesc(imgCapsids_Gauss);axis equal;colormap(gray);

    imgCapsidsMonoInfected=imgCapsids_Gauss;%imfilter(imgCapsids,hG,'replicate');
    figure(210);clf;imagesc(imgCapsidsMonoInfected);axis equal;colormap(gray);

    valIntCaps=imgCapsidsMonoInfected(mskCellMonoInfected>0);
    thld_IntCaps=mean(valIntCaps)+2*std(valIntCaps);
    mskCaps=imgCapsidsMonoInfected>thld_IntCaps;
    imgCapsidsMonoInfected(mskCellMonoInfected==0)=1;
    imgCapsidsMonoInfected(mskCellMonoInfected>0)=2;
    imgCapsidsMonoInfected(mskCellMonoInfected>0&mskCaps>0)=3;
    
    figure(210);clf;imagesc(imgCapsidsMonoInfected);axis equal;colormap(gray);
    figure(21);clf;imagesc(imgCapsids);axis equal;colormap(gray);
    
    % segment using k-means and defined mask of capdsids area as the class
    % with higher intensity
    %nK=3;
    mskLocCapsids=kmeans2D(imgCapsidsMonoInfected,nK_caps);close(50);
    mskLocCapsids=mskLocCapsids==nK_caps;
    %figure(21);clf;imagesc(mskLocCapsids);axis equal;colormap(gray);
    %mskLocCapsids=imclose(mskLocCapsids,strel('disk',1));
    mskLocCapsids=imfill(mskLocCapsids,'holes');
    figure(21);clf;imagesc(mskLocCapsids);axis equal;colormap(gray);
        
%     for iROI=1:nROI
%         if (cellSummary(iROI,3)==1)
%             disp([iROI,cellSummary(iROI,3)])
%             hold on; plot(cellSummary(iROI,6),cellSummary(iROI,7),'m.');
%         end %if
%     end
    
    %% reject small capsids areas
    mskLocCapsids=bwlabel(mskLocCapsids);
    propCapsids=regionprops(mskLocCapsids>0,'Area','Centroid','PixelList');
    nLocCapsids=numel(propCapsids);
    
    areaCaps=zeros(nLocCapsids,1);
    for iL=1:nLocCapsids;areaCaps(iL)=propCapsids(iL).Area;end%for        
    
    removeCapsid=zeros(size(mskLocCapsids));
    for iL=1:nLocCapsids
        if (areaCaps(iL)<minCapsArea)
            removeCapsid(mskLocCapsids==iL)=1;
        end
    end
    mskLocCapsids(removeCapsid>0)=0; 
    mskLocCapsids=mskLocCapsids>0;    
    figure(205);clf;imagesc(mskLocCapsids);axis equal;colormap(gray);
    
    %% reject capsids areas with low intense
    mskLocCapsids=bwlabel(mskLocCapsids);
    propCapsids=regionprops(mskLocCapsids>0,'Area','Centroid','PixelList');
    nLocCapsids=numel(propCapsids);
    
    z=zeros(nLocCapsids,1);
    for iL=1:nLocCapsids
        curCapsidInt=imgCapsids(mskLocCapsids==iL);    
        ratioPixHighInt=mean(curCapsidInt>intROI_neg);
        z(iL)=ratioPixHighInt;            
    end        
    %unique(mskLocCapsids(:))
    for iL=1:nLocCapsids
        if (z(iL)<0.05)
            mskLocCapsids(mskLocCapsids==iL)=0;
        end%if
    end
    %unique(mskLocCapsids(:))
    
    mskLocCapsids=mskLocCapsids>0;
    figure(205);clf;imagesc(mskLocCapsids);axis equal;colormap(gray);
    
    %% Overlay DNA replication and Capsid areas    
    imgDNA=imgInfectionGlobal;%figure(201);clf;imagesc(imgDNA);colormap(gray);% 1 = Bgd; 2= Cells; 3 = DNA Replication
    mskCellInf=1+mskCellMonoInfected; % 1=Bgd + cells with 0 or multiple infections; 2 = cells mono-infected
    mskCellBorder=edge(mskCell>0);
    imgDNA(mskCellInf~=2)=1;% Keep only mono-infected cells    
    mskCellInf(mskLocCapsids>0)=3; % 3 = Capsids area    
    %figure(20);clf;imagesc(mskCellInf);colormap(gray);    
    
    mskOverlap=(mskCellInf==3)&(imgDNA==3);
    imgDNA(mskCellInf==1)=1;
    imgDNA(mskCellInf==3)=4;
    imgDNA(mskOverlap)=5;
    imgDNA(mskCellBorder)=2;
    lstCat=unique(imgDNA(:));
    cmapSeg=[0,0,0;0.5,0.5,0.5;1,0,1;0,1,0;1,1,1];% 1 =bgd ; 2 = cells; 3 = DNA; 4 = capsids; 5 overlay 3 & 5    
    cmapSeg=cmapSeg(lstCat,:); % if there is a missing catergory;
    figure(201);clf;imagesc(imgDNA);axis equal;axis off;colormap(cmapSeg);    
    
    %% Get specifications (localization, shape, ...) of capsids assemblies
    [locCapsids]=get_capsidSpec(mskLocCapsids,mskCell,cellSummary,imgYFP,imgInfectionGlobal,imgH,imgW,nROI);
    
    
% % % % % %     mskLocCapsids=bwlabel(mskLocCapsids>0);
% % % % % %     propCapsids=regionprops(mskLocCapsids>0,'Area','Centroid','PixelList','MajorAxisLength','MinorAxisLength','Perimeter');    
% % % % % %     nLocCapsids=numel(propCapsids);
% % % % % %     
% % % % % %     locCapsids=zeros(nLocCapsids,16);
% % % % % %     shapeCapsids=zeros(nLocCapsids,5);
% % % % % %     
% % % % % %     for iL=1:nLocCapsids
% % % % % %         %shapeCapsids(iL,:)=[bwarea(mskLocCapsids==iL),propCapsids(iL).Perimeter,-1,propCapsids(iL).MajorAxisLength,propCapsids(iL).MinorAxisLength];
% % % % % %         shapeCapsids(iL,:)=[propCapsids(iL).Area,propCapsids(iL).Perimeter,-1,propCapsids(iL).MajorAxisLength,propCapsids(iL).MinorAxisLength];
% % % % % %     end%for
% % % % % %     shapeCapsids(:,3)=shapeCapsids(:,5)./shapeCapsids(:,4);
% % % % % %     
% % % % % %     for iL=1:nLocCapsids
% % % % % %         posCaps=round(propCapsids(iL).Centroid);                
% % % % % %         msk_rows=posCaps(2)-3:posCaps(2)+3;msk_rows(msk_rows<1)=[];msk_rows(msk_rows>size(mskCell,1))=[];
% % % % % %         msk_cols=posCaps(1)-3:posCaps(1)+3;msk_cols(msk_cols<1)=[];msk_cols(msk_cols>size(mskCell,2))=[];
% % % % % %         curROI=mskCell(msk_rows, msk_cols);curROI=curROI(:);curROI(curROI==0)=[];  
% % % % % %         curROI=median(curROI);
% % % % % %         if (mod(curROI,1)~=0) 
% % % % % %             curROI=[];
% % % % % %         end
% % % % % %         
% % % % % %         if ((numel(curROI==1)) & (cellSummary(curROI,3)==1))
% % % % % %             
% % % % % %             locCapsids(iL,1:2)=[iL,curROI];
% % % % % %             
% % % % % %             % mass center of capsid and mean intensity (cols 3-5)
% % % % % %             [pixROI]=propCapsids(iL).PixelList;
% % % % % %             nPix=size(pixROI,1);
% % % % % %             for iPix=1:nPix; pixROI(iPix,3)=imgYFP(pixROI(iPix,2),pixROI(iPix,1));end
% % % % % %             xCM=sum(pixROI(:,1).*pixROI(:,3))/sum(pixROI(:,3));
% % % % % %             yCM=sum(pixROI(:,2).*pixROI(:,3))/sum(pixROI(:,3));
% % % % % %             locCapsids(iL,3:4)=[xCM,yCM];
% % % % % %             locCapsids(iL,5)=mean(pixROI(:,3));
% % % % % %             
% % % % % %             % distance dna replication - capsid assembly (col 6)
% % % % % %             [xyCMdna]=cellSummary(cellSummary(:,1)==locCapsids(iL,2),6:7);
% % % % % %             locCapsids(iL,6)=sqrt((xCM-xyCMdna(1))^2+(yCM-xyCMdna(2))^2);
% % % % % %             
% % % % % %             % overlap btw dna & capsids
% % % % % %             iCell=locCapsids(iL,2);
% % % % % % 
% % % % % %             curCell_DNA=(mskCell==iCell) & (imgInfectionGlobal ==3);%figure(205);clf;imagesc(curCell_DNA);axis equal;colormap(gray);
% % % % % %             curCell_Capsid=mskLocCapsids==iL;%figure(205);clf;imagesc(curCell_Capsid);axis equal;colormap(gray);
% % % % % %             curCell_Overlay=curCell_Capsid & curCell_DNA;%figure(205);clf;imagesc(curCell_Overlay);axis equal;colormap(gray);
% % % % % %             
% % % % % %             [iiDNA,~]=find(curCell_DNA>0);
% % % % % %             [iiCaps,~]=find(curCell_Capsid>0);
% % % % % %             [iiOverlay,~]=find(curCell_Overlay>0);
% % % % % %             
% % % % % %             sDNA=size(iiDNA,1);
% % % % % %             sCaps=size(iiCaps,1);
% % % % % %             sOverlay=size(iiOverlay,1);
% % % % % %             if (sOverlay==0)
% % % % % %                 locCapsids(iL,7)=0;
% % % % % %                 locCapsids(iL,8)=0;
% % % % % %             else
% % % % % %                 locCapsids(iL,7)=100*sOverlay/(sCaps);
% % % % % %                 locCapsids(iL,8)=100*sOverlay/(sDNA);
% % % % % %             end
% % % % % % %            kColor=[0 0 0]/255;
% % % % % % %            bColor=[0 130 254]/255;
% % % % % %             oColor=[255 128 1]/255;
% % % % % % %            mColor=[255 0 255]/255;
% % % % % %             if (locCapsids(iL,7)>0)
% % % % % %                 figure(201);hold on; text(locCapsids(iL,3)-0,locCapsids(iL,4)-10,num2str(locCapsids(iL,7),3),'Color',oColor,'FontName','Arial','FontWeight','bold');
% % % % % %                 figure(201);hold on; plot(locCapsids(iL,3),locCapsids(iL,4),'k+');
% % % % % %             else
% % % % % %                 figure(201);hold on; text(locCapsids(iL,3)+0,locCapsids(iL,4)+10,num2str(locCapsids(iL,7),3),'Color','red','FontName','Arial','FontWeight','bold');
% % % % % %                 figure(201);hold on; plot(locCapsids(iL,3),locCapsids(iL,4),'ko');
% % % % % %             end          
% % % % % % 
% % % % % %             %%
% % % % % %             j0=xCM;j1=xyCMdna(1);i0=yCM;i1=xyCMdna(2);
% % % % % %             imDNA=bwlabel(imgInfectionGlobal==3);
% % % % % %             curIDdna=imDNA(round(xyCMdna(2)),round(xyCMdna(1)));            
% % % % % %             
% % % % % %             imProfile=zeros(size(mskCell));
% % % % % %             imProfile(imDNA==curIDdna)=1;
% % % % % %             imProfile(mskLocCapsids==iL)=2;            
% % % % % %             
% % % % % % %             figure(21);clf;imagesc(imProfile);axis equal;colormap(gray);
% % % % % % %             figure(21);hold on; plot([j0,j1],[i0,i1],'y-.');
% % % % % %             jProfile=j0:sign(j1-j0):j1;
% % % % % %             iProfile=i0:sign(i1-i0):i1;
% % % % % %             
% % % % % %             if(numel(iProfile)>numel(jProfile))
% % % % % %                 jProfile=linspace(jProfile(1),jProfile(end),numel(iProfile));
% % % % % %             else
% % % % % %                 iProfile=linspace(iProfile(1),iProfile(end),numel(jProfile));
% % % % % %             end
% % % % % %             jProfile=floor(0.5+jProfile);
% % % % % %             iProfile=floor(0.5+iProfile);
% % % % % %             
% % % % % %             jProfile((iProfile>imgH) | (iProfile<1))=[];
% % % % % %             iProfile((iProfile>imgH) | (iProfile<1))=[];
% % % % % %             iProfile((jProfile>imgW) | (jProfile<1))=[];
% % % % % %             jProfile((jProfile>imgW) | (jProfile<1))=[];
% % % % % %             %hold on; plot(jProfile,iProfile,'w.');
% % % % % %             lProfile=sqrt((iProfile(1)-iProfile(end))^2 + (jProfile(1)-jProfile(end))^2);
% % % % % %             nP=numel(iProfile);
% % % % % %             intCell=zeros(nP,3);
% % % % % %             intCell(:,1)=linspace(0,lProfile,nP);
% % % % % %             intCell(:,3)=1:nP;
% % % % % %             %intGP=intSPP1;
% % % % % %             for iP=1:nP
% % % % % %                 intCell(iP,2)=imProfile(iProfile(iP),jProfile(iP));
% % % % % %                 %   intGP(iP,2)=imgYFP_profile(iProfile(iP),jProfile(iP));
% % % % % %             end%for
% % % % % %             intCell(intCell(:,2)~=0,:)=[];
% % % % % %             
% % % % % %             if (~isempty(intCell))
% % % % % %                 locCapsids(iL,9)=intCell(end,1)-intCell(1,1);
% % % % % %             else
% % % % % %                 locCapsids(iL,9)=0;
% % % % % %             end
% % % % % %             %%
% % % % % %             
% % % % % %             %
% % % % % %             xCaps=locCapsids(iL,3);
% % % % % %             yCaps=locCapsids(iL,4);            
% % % % % %             xRep=cellSummary(locCapsids(iL,2),6);
% % % % % %             yRep=cellSummary(locCapsids(iL,2),7);            
% % % % % %             xPole=cellSummary(locCapsids(iL,2),14);
% % % % % %             yPole=cellSummary(locCapsids(iL,2),15);
% % % % % %             
% % % % % %             
% % % % % %             xPole=xPole-xRep;
% % % % % %             xCaps=xCaps-xRep;
% % % % % %             if (xPole<0)
% % % % % %                 xPole=-xPole;
% % % % % %                 xCaps=-xCaps;
% % % % % %             end
% % % % % %             
% % % % % %             yPole=yPole-yRep;
% % % % % %             yCaps=yCaps-yRep;
% % % % % %             if (yPole<0)
% % % % % %                 yPole=-yPole;
% % % % % %                 yCaps=-yCaps;
% % % % % %             end
% % % % % %             
% % % % % %             %figure(900);clf;hold on;plot(0,0,'go');plot(xPole,yPole,'ks');plot(xCaps,yCaps,'r+');axis equal            
% % % % % %             
% % % % % %             if (abs(xCaps)>abs(yCaps))
% % % % % %                 testPos=xCaps;
% % % % % %             else
% % % % % %                 testPos=yCaps;
% % % % % %             end
% % % % % %             if testPos>0;
% % % % % %                 locCapsids(iL,16)=1;
% % % % % %             else
% % % % % %                 locCapsids(iL,16)=2;
% % % % % %             end
% % % % % %             %title(num2str(locCapsids(iL,16)));pause()
% % % % % %             
% % % % % %         else
% % % % % %             
% % % % % %             locCapsids(iL,2:end)=-1;
% % % % % %         end
% % % % % %     end%for
% % % % % %     
% % % % % %     for iROI=1:nROI
% % % % % %         distZ=locCapsids(locCapsids(:,2)==iROI,6);
% % % % % %         [z,indS]=sort(distZ);
% % % % % %         %[(1:numel(indS))',indS]
% % % % % %         nS=numel(indS);inv_indS=zeros(nS,1);
% % % % % %         for iS=1:nS;inv_indS(indS(iS))=iS;end
% % % % % %         locCapsids(locCapsids(:,2)==iROI,10)=inv_indS;
% % % % % %     end
% % % % % %     
% % % % % %     locCapsids(:,11:15)=shapeCapsids;
    
end%if
locCapsids(locCapsids(:,1)==0,:)=[];
disp('Capsid loc parameters done');
end%function

