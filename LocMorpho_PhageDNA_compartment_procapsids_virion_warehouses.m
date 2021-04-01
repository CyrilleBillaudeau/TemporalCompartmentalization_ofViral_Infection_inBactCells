function []=LocMorpho_PhageDNA_compartment_procapsids_virion_warehouses()
%% ========================================================================
%% This scripts aims to detect and indentify areas associated to SPP1 phage
%% DNA, viral encapsidation (pro-capsids) and virion warehouse (capsids).
%% During step #1 '*_mask0.txt' files has been created, and has to be
%% checked, corrected (if needed) and saved in a new file named '*_mask.txt'.
%%
%% Here the script will load the file list created in step #1 to load
%% raw data and binary mask and perform an automatic detection of areas
%% associated to phage DNA replication and the capsid formations.
%% Results are saved in several files:
%% **** for each file of the list
%% - parameters used in analysis are saved in a file named as the raw data
%% by replacing '.tif' with '_paramSegmentation_RepCaps.txt'.
%% - figures displaying final segmentation and quantification on SPP1 phage
%% DNA replication and (pro-)capsids areas are saved in folder 'Fig' with
%% respective names as the raw data by replacing '.tif' with
%% '_seg.tif' (or '_seg.png)', '_q1_SPP1_PhageDNA_LocMorphology.png' and
%% '_q2_ProCapsid_LocMorphology.png'.
%% - a text image in 'Mask' folder showing binary mask of cell and Phage DNA
%% **** for the compilation of all files of the list in the 'Resume' folder:
%% - all quantification performed on the file list are saved in two
%% tables named as the filename used for the list by replacing '.mat' with
%% '_tabSummary.txt' and '_tabSummaryCapsids.txt'
%% - figures representing some quantifications are saved in in two figures
%% named as the filename used for the list by replacing '.mat' with
%% '_stat.png' and '_statCapsids.png'
%%
%% Required: Image Processing Toolbox
%%
%% -- Cyrille Billaudeau / INRAE / Micalis / ProCeD Lab / 210209
%% ------------------------------------------------------------------------

close all
clear all

%% Initialization
% add matlab functions to Matlab Path (alternatively, you can directly set
% your path and comment the next four lines).
scriptPath=mfilename('fullpath');
scriptPath=scriptPath(1:max(strfind(scriptPath,filesep)));
rootPathFunctions=strcat([scriptPath,'functions']);
addpath(genpath(rootPathFunctions));

% path contaning data: next lines has to be modified according to your data
% organization (location). 
cd(scriptPath);
cd ../Data
%startPath='E:\path\to\your\Data';cd(startPath);

% select file liste created at step #1
[fileLst,pathLst]= uigetfile('select file list');

% Intialization of parameters used in pro-capsids and capsids localization
doGlobalSegForSPP1phageDNA=0;% Define if SPP1 phade DNA segmentation should be done by taking account all files or by single files
limAreaInf=7;% Min area for DNA replication area (pix unit)
nK_caps=3;% Number of cluster used in capsid detection
minCapsArea=6; % Min area of capsid detection (pix unit)
isCapsid=1; % Define if detection concerns pro-capsid (0) or capsid (1)

prompt = {...
    'SPP1 phage DNA segmentation: single (0) / global (1)',...
    'Minimal area for SPP1 phage DNA replication (pixel unit):',...
    'Number of cluster groups in capsid detection:',...
    'Minimal area for capsid (pixel unit):',...
    'Detection of pro-capsids (0) / capsids (1)'};
dlg_title = 'Parameters';
num_lines = 1;
defaultans = {num2str(doGlobalSegForSPP1phageDNA),num2str(limAreaInf),num2str(nK_caps),num2str(minCapsArea),num2str(isCapsid)};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
doGlobalSegForSPP1phageDNA=str2num(answer{1});
limAreaInf=str2num(answer{2});
nK_caps=str2num(answer{3});
minCapsArea=str2num(answer{4});
isCapsid=str2num(answer{5});

cd(pathLst);lstFolder=load(fileLst);lstFolder=lstFolder.lstFolder;
nFolder=numel(lstFolder);

%% Find all SPP1 DNA areas in seperate folder based on kmeans segmentation (nK=3)
%% and extract fluorescence intensity (min & max) of areas associated to
%% background, cytoplasm (outside or wo/ SPP1 DNA) and SPP1 DNA
if (doGlobalSegForSPP1phageDNA==1)
    
    [tab_ThldInf]=estimateActivationThreshold(lstFolder,nFolder);
    
    % Indentify files were SPP1 DNA segmentation was detected to get correct intensity threshold.
    doFindOutlier=1;
    if (doFindOutlier)
        tab_ThldInf(tab_ThldInf(:,1)==0,:)=NaN;
        figure(980);clf;hold on;xlim([0 nFolder+1])
        plot(tab_ThldInf(:,1),tab_ThldInf(:,2),'--o');
        plot(tab_ThldInf(:,1),tab_ThldInf(:,3),'--o');
        plot(tab_ThldInf(:,1),tab_ThldInf(:,4),'--o');
        plot(tab_ThldInf(:,1),tab_ThldInf(:,5),'--o');
        plot(tab_ThldInf(:,1),tab_ThldInf(:,6),'--o');
        plot(tab_ThldInf(:,1),tab_ThldInf(:,7),'--o');
        xlabel('Folder #');ylabel('Fluo intensity (a.u.)');
        legend('min BGD','max BGD','min cytoplasm (outside or wo/ SPP1 DNA)','max  (outside or wo/ SPP1 DNA)','min SPP1 DNA','max  SPP1 DNA','Location','northoutside')
        %valThld=median(tab_ThldInf(:,6));
        
        outlier1=(tab_ThldInf(:,6)>nanmean(tab_ThldInf(:,6))+1.5*nanstd(tab_ThldInf(:,6)));% find folder presenting outlier max value for cell (oustide or wo/ SPP1 DNA)
        outlier2=(tab_ThldInf(:,7)>nanmean(tab_ThldInf(:,7))+1.5*nanstd(tab_ThldInf(:,7)));% find folder presenting outlier min value for SPP1 DNA area
        outlier=outlier1&outlier2
        tab_ThldInf(outlier,6)=nanmedian(tab_ThldInf(~outlier,6)); % change threshold of outlier files by median of non-outlier file
        if sum(outlier>0)
            plot(tab_ThldInf(:,1),tab_ThldInf(:,6),'k-.o');
        end%if
    end%if(doFindOutlier)
end%if

%% Detection, localization and morphology of areas associated to the
%% replication of SPP1 phage DNA, and the assembly and warehouse of
%% virion (respectively pro-caspid and capsid)

tab_cellSummary=[];
tab_locCapsids=[];

hWB = waitbar(0,'File are processed ...');
for iFolder=1:nFolder
    disp('===============================================================')
    % Get filename for raw acquisitions and binary mask
    path=lstFolder{iFolder};
    indexSep=1+max(strfind(lstFolder{iFolder},filesep));
    imgFilename=path(indexSep:end);
    pathImg=path(1:(indexSep-2));
    disp(pathImg)
    disp(imgFilename)
    fileGen=imgFilename(1:strfind(imgFilename,'.tif')-1);
    cd(pathImg)
    cd('Output/');
    mskFilename=strcat([fileGen,'_mask.txt']);
    
    if (exist(mskFilename, 'file')==0)
        disp(strcat(['/!\ mask not found for: ',mskFilename]))
        cd ../
    else
        %% Load Images and data pre-processing
        mskCell=load(strcat([fileGen,'_mask.txt']));
        mskCell=mskCell>0;
        mskCell=bwlabel(mskCell);
        cd ../
        
        imgInfo=imfinfo(imgFilename);
        nChannel=numel(imgInfo);
        
        imgW=imgInfo(1).Width;
        imgH=imgInfo(1).Height;
        pixelSize=1/imgInfo(1).XResolution;% µm 0.1024
        
        imgYFP=double(imread(imgFilename,1)); % GP
        imgRFP=double(imread(imgFilename,2)); % SPP1 DNA
        imgBF=double(imread(imgFilename,3));
        
        imgYFP_profile=imgYFP;imgYFP_profile(mskCell==0)=NaN;
        imgRFP_profile=imgRFP;imgRFP_profile(mskCell==0)=NaN;
        
        propROI=regionprops(mskCell,'Area','Centroid','Orientation','BoundingBox','PixelList');
        numROI=numel(propROI);
        figure(16);clf;imagesc(mskCell);cmap=[0,0,0;lines(numROI)];colormap(cmap);
        valSPP1=[];
        valSPP1_cell=[];
        nROI=numel(propROI);
        cellSummary=zeros(nROI,15);
        for iROI=1:nROI
            cellSummary(iROI,1)=iROI;
            cellSummary(iROI,2)=propROI(iROI).Area;
        end%for
        
        %% Localization and morphology of the phage DNA compartment
        if (doGlobalSegForSPP1phageDNA~=1)
            [tabInfectedCells,imgInfectionGlobal,imgInfectionGlobalHigh,propROI_infection]=findInfectedCells_single(imgRFP,mskCell,limAreaInf);
        else
            valThld=tab_ThldInf(iFolder,6)
            [tabInfectedCells,imgInfectionGlobal,imgInfectionGlobalHigh,propROI_infection]=findInfectedCells_global(imgRFP,mskCell,limAreaInf,valThld);
        end
        
        cellSummary(:,3)=tabInfectedCells;
        cellSummary=extractTotInt_InfectedCells(imgRFP,mskCell,cellSummary);
        cellSummary=extractParam_DNA_replication(cellSummary,propROI,imgRFP,imgInfectionGlobal);
        cellSummary=extractLoc_DNA_replication(cellSummary,mskCell,propROI,imgH,imgW,pixelSize);
        
        %% Localization and morphology of procapsids and virion warehouses
        if (isCapsid)
            [locCapsids]=extractLocCapsids(cellSummary,mskCell,imgYFP,imgInfectionGlobal,minCapsArea,imgW,imgH,imgFilename);
        else
            [locCapsids]=extractLocProCapsids(cellSummary,mskCell,imgYFP,imgInfectionGlobal,minCapsArea,imgW,imgH,nK_caps);
        end%if
        %pause();
        
        %% Export graphics representing:
        %% - localization and morphology of areas detected for the
        %% replication of SPP1 phage DNA, viral encapsidation (pro-capsid)
        %% and virion warehouse (capsid). Two image format (png and tiff)
        %% are proposed.
        %% - localization and morphology of SPP1 phage DNA
        %% - Quantifications on SPP1 phage DNA
        %% - Quantifications on (pro-)capsid
        
        cd(pathImg)
        cd('Output');
        fileGen=imgFilename(1:strfind(imgFilename,'.tif')-1);
        paramSeg=[doGlobalSegForSPP1phageDNA;limAreaInf;nK_caps;minCapsArea;isCapsid];
        save(strcat([fileGen,'_paramSegmentation_RepCaps.txt']),'paramSeg','-ascii');
        
        % replication, pro-capsids and capsid
        exportImageInTiff=1;
        cd Fig
        fig=gcf;fig.PaperOrientation='portrait';fig.PaperType='a4';fig.Units='Normalized';fig.OuterPosition=[0 0 1 1];axis equal;
        if (exportImageInTiff)
            imgOut=strcat([imgFilename(1:end-4),'_seg.tif']);
            saveas(201,imgOut,'tiff');
        else
            imgOut=strcat([imgFilename(1:end-4),'_seg.png']);
            print(201,imgOut,'-dpng','-r300');
        end
        cd ../../
        
        % replication areas
        figure(18);clf;imagesc(imgInfectionGlobal);colormap(gray);
        for iROI=1:nROI
            posROI=propROI(iROI).Centroid;
            text(posROI(1)-20,posROI(2)-20,num2str(iROI),'Color','w')
            hold on; plot(cellSummary(iROI,6),cellSummary(iROI,7),'m+');
        end
        
        if (~exist('Output/Mask','dir'))
            mkdir('Output/Mask');
        end%if
        
        cd Output/Mask/
        imgOut=strcat([imgFilename(1:end-4),'_mask_DNA_Cell.txt']);
        save(imgOut,'imgInfectionGlobal','-ascii');
        cd ../../
        
        
        if (~isempty(cellSummary))
            
            tab_cellSummary=[tab_cellSummary;iFolder*ones(size(cellSummary,1),1),cellSummary];
            
            f1=figure(100);clf;hold on;
            f1.Units='normalized';
            f1.Position=[0 0 1. 1.];
            szFtTitle=10;
            szFtAxis=szFtTitle;
            
            [n,xout]=hist(cellSummary(:,3),[0:3]);n=100*n/sum(n);
            subplot(3,2,1);bar(xout,n,'k');xlim([-0.5 3.5]);ylim([0 100]);xlabel('SPP1 Phage DNA / cell','FontSize',szFtAxis);title('Distribution of SPP1 Phg DNA per cell (%)','FontSize',szFtTitle)
            
            cellSummary=cellSummary(cellSummary(:,3)==1,:); % remove non-mono infected cells
            [n,xout]=hist(cellSummary(:,8)*pixelSize^2,[0:0.1:2.5]);n=100*n/sum(n);
            subplot(3,2,2);bar(xout,n,'k');xlim([-xout(2) xout(end)+xout(2)]);xlabel('Cell area (um2)','FontSize',szFtAxis);title('Area of SPP1 Phg DNA','FontSize',szFtTitle)
            
            [n,xout]=hist(cellSummary(:,9)*pixelSize,[0:0.1:2.5]);n=100*n/sum(n);
            subplot(3,2,3);bar(xout,n,'k');xlim([-xout(2) xout(end)+xout(2)]);xlabel('Length of major axis (um)','FontSize',szFtAxis);title('Length of SPP1 Phg DNA in mono-infected cells','FontSize',szFtTitle)
            
            [n,xout]=hist(cellSummary(:,11),[0:20:2000]);n=100*n/sum(n);
            subplot(3,2,4);bar(xout,n,'k');xlim([-xout(2) xout(end)+xout(2)]);xlabel('Fluorescence intensity (a.u.)','FontSize',szFtAxis);title('Int. of SPP1 Phg DNA in mono-infected cells','FontSize',szFtTitle)
            
            [n,xout]=hist(cellSummary(:,12)*pixelSize,[0:0.2:4]);n=100*n/sum(n);
            subplot(3,2,5);bar(xout,n,'k');xlim([-xout(2) xout(end)+xout(2)]);xlabel('Distance to the proximal pole (µm)','FontSize',szFtAxis);;title('Dist. of the SPP1 Phg DNA to the proximal cell pole in mono-infected cells','FontSize',szFtTitle);
            
            [n,xout]=hist(cellSummary(:,12)./cellSummary(:,13),[0:0.025:0.5]);n=100*n/sum(n);
            subplot(3,2,6);bar(xout,n,'k');xlim([-xout(2) xout(end)+xout(2)]);xlabel('Normalized distance to the proximal pole (µm)','FontSize',szFtAxis);;title('Normalized Dist. of the SPP1 Phg DNA to the proximal cell pole in mono-infected cells','FontSize',szFtTitle);
            
            cd Output/Fig
            imgOut=strcat([imgFilename(1:end-4),'_q1_SPP1_PhageDNA_LocMorphology.png']);
            print(100,imgOut,'-dpng','-r300');
            cd ../../
            
        end
        
        if (~isempty(locCapsids))
            tab_locCapsids=[tab_locCapsids;iFolder*ones(size(locCapsids,1),1),locCapsids];
            maxCapsids=max(locCapsids(:,10));
            figure(101);clf;
            subplot(2,2,1);plot(locCapsids(:,10),locCapsids(:,6)*pixelSize,'k.');xlim([0 maxCapsids+1])
            xlabel('rank with increasing distance');ylabel('distance replication - (pro)capsids (um)')
            subplot(2,2,2);plot(locCapsids(:,10),locCapsids(:,5),'k.');xlim([0 maxCapsids+1])
            xlabel('rank with increasing distance');ylabel('mean intensity (pro)capsids (a.u.)')
            subplot(2,2,3);plot(locCapsids(:,10),locCapsids(:,7),'k.');xlim([0 maxCapsids+1])
            xlabel('rank with increasing distance');ylabel('overlap (pro)capsids / replication (%)')
            subplot(2,2,4);plot(locCapsids(:,10),locCapsids(:,8),'k.');xlim([0 maxCapsids+1])
            xlabel('rank with increasing distance');ylabel('overlap replication - (pro)capsids (%)')
            
            cd Output/Fig
            if (isCapsid)
                imgOut=strcat([imgFilename(1:end-4),'_q2_Capsid_LocMorphology.png']);
            else
                imgOut=strcat([imgFilename(1:end-4),'_q2_ProCapsid_LocMorphology.png']);
            end
            print(101,imgOut,'-dpng','-r300');
            cd ../../
        end
        
        
    end%if msk exist
    waitbar(iFolder / nFolder);
    %pause();
end%for iFolder
close(hWB);
close(201);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Plot general outputs on SPP1 phage DNA & Table creation - Saving %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% SPP1 Phage DNA
% - Distribution of SPP1 Phage DNA replication foci per cell
% - Cell area occupied by SPP1 Phage DNA (um2)
% - Length of major axis of area occupied by SPP1 Phage DNA (um)
% - Average intensity of fluorescence associated to SPP1 Phage DNA
% - Distance of the SPP1 DNA focus center of mass to the proximal cell pole
% - Normalized distance of the SPP1 DNA focus center of mass to the
% proximal cell pole to total cell length.

cellSummary=tab_cellSummary(:,2:end);
f1=figure(100);clf;hold on;
f1.Units='normalized';
f1.Position=[0 0 1. 1.];
szFtTitle=10;
szFtAxis=szFtTitle;

[n,xout]=hist(cellSummary(:,3),[0:3]);n=100*n/sum(n);
subplot(3,2,1);bar(xout,n,'k');xlim([-0.5 3.5]);ylim([0 100]);xlabel('SPP1 Phage DNA / cell','FontSize',szFtAxis);title('Distribution of SPP1 Phg DNA per cell (%)','FontSize',szFtTitle)

cellSummary=cellSummary(cellSummary(:,3)==1,:); % remove non-mono infected cells
[n,xout]=hist(cellSummary(:,8)*pixelSize^2,[0:0.1:2.5]);n=100*n/sum(n);
subplot(3,2,2);bar(xout,n,'k');xlim([-xout(2) xout(end)+xout(2)]);xlabel('Cell area (um2)','FontSize',szFtAxis);title('Area of SPP1 Phg DNA','FontSize',szFtTitle)

[n,xout]=hist(cellSummary(:,9)*pixelSize,[0:0.1:2.5]);n=100*n/sum(n);
subplot(3,2,3);bar(xout,n,'k');xlim([-xout(2) xout(end)+xout(2)]);xlabel('Length of major axis (um)','FontSize',szFtAxis);title('Length of SPP1 Phg DNA in mono-infected cells','FontSize',szFtTitle)

[n,xout]=hist(cellSummary(:,11),[0:50:2000]);n=100*n/sum(n);
subplot(3,2,4);bar(xout,n,'k');xlim([-xout(2) xout(end)+xout(2)]);xlabel('Fluorescence intensity (a.u.)','FontSize',szFtAxis);title('Int. of SPP1 Phg DNA in mono-infected cells','FontSize',szFtTitle)

[n,xout]=hist(cellSummary(:,12)*pixelSize,[0:0.2:4]);n=100*n/sum(n);
subplot(3,2,5);bar(xout,n,'k');xlim([-xout(2) xout(end)+xout(2)]);xlabel('Distance to the proximal pole (µm)','FontSize',szFtAxis);;title('Dist. of the SPP1 Phg DNA to the proximal cell pole in mono-infected cells','FontSize',szFtTitle);

[n,xout]=hist(cellSummary(:,12)./cellSummary(:,13),[0:0.025:0.5]);n=100*n/sum(n);
subplot(3,2,6);bar(xout,n,'k');xlim([-xout(2) xout(end)+xout(2)]);xlabel('Normalized distance to the proximal pole (µm)','FontSize',szFtAxis);;title('Normalized Dist. of the SPP1 Phg DNA to the proximal cell pole in mono-infected cells','FontSize',szFtTitle);

%% Pro-capsids & Virion warehouse (capsids)
% - xAxis: detected (pro-)capsids in mono-infected cells are sorted
% according to the distance between SPP1 phage DNA area and capsid foci
% (center to center). The resulting rank is then provided from nearest
% ('1') to farest (>1).
% - boxplot corresponding to overlapping areas (between SPP1 phage DNA and
% (pro-)capsid: non-overlapping area (i.e. there is a distance gap between
% those two areas) are not taking into account.

locCapsids=tab_locCapsids(:,2:end);
if (~isempty(locCapsids))
    
    f2=figure(101);clf;
    f2.Units='normalized';
    f2.Position=[0 0 1. 1.];
    szFtTitle=10;
    szFtAxis=szFtTitle;
    
    maxCapsids=max(locCapsids(:,10));
    sizeRank=0;
    for iM=1:maxCapsids;sizeRank=max([sizeRank;sum(locCapsids(:,10)==iM)]);end
    
    dataPlot=NaN(sizeRank,maxCapsids);
    for iM=1:maxCapsids;
        ind=locCapsids(:,10)==iM;
        dataPlot(1:sum(ind),iM)=locCapsids(ind,6)*pixelSize;
    end
    subplot(2,2,1);boxplot(dataPlot,'colors','k','symbol','k.');xlim([0 maxCapsids+1])
    xlabel('Foci rank according to their distance from the replication area center','FontSize',szFtAxis);
    ylabel('distance replication - (pro)capsids (um)','FontSize',szFtAxis);
    title('Distance in mono-infected cells','FontSize',szFtTitle);
    
    dataPlot=NaN(sizeRank,maxCapsids);
    for iM=1:maxCapsids;ind=locCapsids(:,10)==iM;dataPlot(1:sum(ind),iM)=locCapsids(ind,5);end
    subplot(2,2,2);boxplot(dataPlot,'colors','k','symbol','k.');xlim([0 maxCapsids+1])
    xlabel('Foci rank according to their distance from the replication area center','FontSize',szFtAxis);
    ylabel('Average intensity (pro)capsids (a.u.)','FontSize',szFtAxis);
    title('Capsids intensity in mono-infected cells','FontSize',szFtTitle);
    
    dataPlot=NaN(sizeRank,maxCapsids);
    for iM=1:maxCapsids;ind=locCapsids(:,10)==iM;dataPlot(1:sum(ind),iM)=locCapsids(ind,7);dataPlot(dataPlot(:,iM)==0,iM)=NaN;end
    subplot(2,2,3);boxplot(dataPlot,'colors','k','symbol','k.');xlim([0 maxCapsids+1]);ylim([0 100])
    xlabel('Foci rank according to their distance from the replication area center','FontSize',szFtAxis);
    ylabel('Overlap (pro)capsids / replication (%)','FontSize',szFtAxis);
    title('Areas overlap in mono-infected cells and with overlap','FontSize',szFtTitle);
    
    dataPlot=NaN(sizeRank,maxCapsids);
    for iM=1:maxCapsids;ind=locCapsids(:,10)==iM;dataPlot(1:sum(ind),iM)=locCapsids(ind,8);dataPlot(dataPlot(:,iM)==0,iM)=NaN;end
    subplot(2,2,4);boxplot(dataPlot,'colors','k','symbol','k.');xlim([0 maxCapsids+1]);ylim([0 100])
    xlabel('Foci rank according to their distance from the replication area center','FontSize',szFtAxis);
    ylabel('Overlap replication - (pro)capsids (%)','FontSize',szFtAxis);
    title('Areas overlap in mono-infected cells and with overlap','FontSize',szFtTitle);
    
end%if

%% Save tables and created graphics in Resume folder:

cd(startPath);
if (~exist('Resume','dir'))
    mkdir('Resume');
end
cd('Resume');

fileOut=fileLst(1:end-4)
save(strcat([fileOut,'_tabSummary.txt']),'tab_cellSummary','-ascii');
print(100,strcat([fileOut,'_stat.png']),'-dpng');
if (~isempty(locCapsids))
    fileOut=fileLst(1:end-4)
    save(strcat([fileOut,'_tabSummaryCapsids.txt']),'tab_locCapsids','-ascii');
    print(101,strcat([fileOut,'_statCapsids.png']),'-dpng');
end%if

cd(startPath);
disp('Done!')

end%function