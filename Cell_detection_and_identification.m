function []=Cell_detection_and_identification()
%% ========================================================================
%% This scripts aims to detect and indentify bacterial cells and return them
%% a binary mask. This first step is not perfect and need to be checked and
%% manually corrected to separe filamentous cells (no membrane staining
%% available for this acquisition and analysis).
%% Results are saved in a folder 'Output' located in the beside the analyzed
%% images:
%% - An image text file (named as the raw data by replacing '.tif' by
%% '_mask0.txt') with the segmented cells. This file has to be checked and
%% corrected (e.g. with ImageJ/Fiji) and saved with same filename (by replacing
%% 'mask0.txt' with 'mask.txt'
%% - All ROI (cells or chained cells) are save in 'Fig' folder and named as
%% the raw data by replacing '.tif' by '_ROI.png')
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

% Edit list of file to process: create new one and save it (fileLstExist=0)
% or load previous one (1).
fileLstExist=1;
if (fileLstExist==0)
    lstFolder=uipickfiles('FilterSpec','*.tif');
    nFolder=numel(lstFolder);
    
    answer = inputdlg('e.g.: control','set filename for save your filelist');
    lstOutput=strcat([cell2mat(answer),'.mat']);
    save(lstOutput,'lstFolder');
else
    [fileLst,pathLst]= uigetfile('select file list')
    cd(pathLst);lstFolder=load(fileLst);lstFolder=lstFolder.lstFolder;
    nFolder=numel(lstFolder);
end

tabROI_folder=zeros(nFolder,3);
doSizeFiltering=0;

for iFolder=1:nFolder
    
    disp('##############################"#')
    path=lstFolder{iFolder};
    indexSep=1+max(strfind(lstFolder{iFolder},filesep));
    imgFilename=path(indexSep:end);
    pathImg=path(1:(indexSep-2));
    disp(pathImg);
    disp(imgFilename);
    
    %% Load acquisition
    cd(pathImg)
    imgInfo=imfinfo(imgFilename);
    nChannel=numel(imgInfo);
    
    imgW=imgInfo(1).Width;
    imgH=imgInfo(1).Height;
    pixelSize=1/imgInfo(1).XResolution;% µm
    
    imgYFP=double(imread(imgFilename,1)); % GP
    imgRFP=double(imread(imgFilename,2)); % SPP1 DNA
    imgBF=double(imread(imgFilename,3));
    
    %% Cell segmentation
    imgYFP=imtophat(imgYFP,strel('disk',15)); %% correct uneven illumination
    figure(10);clf;imagesc(imgYFP);axis equal;colormap(gray);title('YFP Channel')
    
    % segment cells based on kmeans segmentation. cell separation is not
    % possible without membrane labeling
    nCluster=2;
    mask1=kmeans2D(double(medfilt2(imgYFP, [3 3])),nCluster);close(50);
    cellMask_yfp=(mask1==nCluster);figure(101);clf;imagesc(cellMask_yfp);axis equal;colormap(gray);
    distrib1stKM=zeros(nCluster,1);
    for iCluster=1:nCluster; distrib1stKM(iCluster)=mean(mask1(:)==iCluster);end
    
    % segment background based on YFP images
    limClusterBGD=1;
    mask1=imerode(mask1<=limClusterBGD,strel('disk',3));%
    figure(11);clf;imagesc(mask1);axis equal;colormap(gray);title('BG pixels by kmeans')
    
    % segment cell using BF (after correction) and using segmentation done on
    % YFP to remove BGD. Objective is to separate cell contour from inside and
    % BGD.
    background = imopen(imgBF,strel('disk',15));
    background = medfilt2(background, [15 15]);
    %figure(12);clf; imagesc(background); axis equal; colormap(gray)
    imgBFcorr=imgBF-background;
    figure(12);clf;
    subplot(2,1,1);imagesc(imgBF); axis equal; colormap(gray);title('BF')
    subplot(2,1,2);imagesc(imgBFcorr); axis equal; colormap(gray);title('BF corr=(BF-imopen(BF))')
    
    imDOG1=imfilter(imgBFcorr,fspecial('gaussian', 25, 2.5),'replicate');
    imDOG2=imfilter(imgBFcorr,fspecial('gaussian', 25, 2.6),'replicate');
    imDOG=imDOG1-imDOG2;
    figure(13);clf;imagesc(imDOG);axis equal;colormap(gray);title('DoG(BF corr)')
    % valBF=imDOG(:);%valBF(valBF==0)=[];
    % figure(14);clf;hist(valBF,100);title('intensities of DoG(BF corr)')
    imDOG(mask1)=-100;
    mask2=kmeans2D(imDOG,3);title('kmeans on DoG(BF corr) taking into consideration 1st kmeans')%close(50);
    
    mskCell=mask2==3;
    mskCell=bwlabel(mskCell);
    %propROI=regionprops(mskCell,'Area');
    nROI=max(mskCell(:));
    tabROI_folder(iFolder,1)=nROI;
    disp(strcat(['# ROI founds before rejecting size filtering: ',num2str(nROI)]))
    
    if (doSizeFiltering)
        areaROI_min=150*(pixelSize^2);%default: 250
        areaROI_min_pix=areaROI_min/(pixelSize^2);
        for iROI=1:nROI
            pixelFound=mskCell==iROI;
            if (sum(pixelFound(:))<areaROI_min_pix)
                mskCell(pixelFound)=0;
                %sum(pixelFound(:))
            end
        end%for
        mskCell=bwlabel(mskCell>0);
        mskCell=imdilate(mskCell,strel('disk',1));
        propROI=regionprops(mskCell,'Area','Centroid');
        nROI=numel(propROI);
        tabROI_folder(iFolder,2)=nROI;
        disp(strcat(['# ROI founds after rejecting small ones: ',num2str(nROI)]))
        
        % Reject big ROIs (cell badly segmented or chained).
        figure(150);clf;imagesc(mskCell);cmap=[0,0,0;lines(nROI)];colormap(cmap);
        areaROI_max=7;
        roiArea=zeros(nROI,1);
        for iROI=1:nROI
            roiArea(iROI)=propROI(iROI).Area*(pixelSize^2);
            xy=propROI(iROI).Centroid;
            if roiArea(iROI)>areaROI_max
                text(xy(1)+5,xy(2)+5,num2str(iROI),'Color','w')
            end
        end
        title('numbered ROI are rejected from analysis');
        
        % keep ROI with ROI
        numROI=0;
        mskCell_tmp=zeros(size(mskCell));
        for iROI=1:nROI
            if roiArea(iROI)<areaROI_max
                numROI=numROI+1;
                mskCell_tmp=mskCell_tmp+(mskCell==iROI)*numROI;
            end
        end
        mskCell=mskCell_tmp;
        figure(16);clf;imagesc(mskCell);cmap=[0,0,0;lines(numROI)];colormap(cmap);
        tabROI_folder(iFolder,3)=nROI;
        disp(strcat(['# ROI kept for analysis: ',num2str(numROI)]))
        
        close(15)
        
    else
        figure(16);clf;imagesc(mskCell);cmap=[0,0,0;lines(nROI)];colormap(cmap);
        tabROI_folder(iFolder,3)=nROI;
        disp(strcat(['# ROI kept for analysis (without size filtering): ',num2str(nROI)]))
        areaROI_min=-1;areaROI_max=-1;
    end
    
    % Export segmentation output
    if (~exist('Output','dir'))
        mkdir('Output');
    end
    
    cd('Output');
    fileGen=imgFilename(1:strfind(imgFilename,'.tif')-1);
    paramSeg=[pixelSize;nCluster;limClusterBGD;areaROI_min;areaROI_max;doSizeFiltering];
    save(strcat([fileGen,'_paramSegmentation.txt']),'paramSeg','-ascii');
    if (~exist('Fig','dir'))
        mkdir('Fig');
    end%if
    cd('Fig');
    figure(16);set(gcf,'units','normalized','outerposition',[0 0 imgH/imgW 1]);
    print(16,strcat([fileGen,'_ROI.png']),'-dpng');
    cd ../
    save(strcat([fileGen,'_mask0.txt']),'mskCell','-ascii');
    disp('##### End of segmentation ######')
    
    close(10);
    close(11);
    close(12);
    close(13);
    close(50);
    %close(15);
    close(101);
    
end

disp('##### End of step #1 ######');
%% --------------------------------------
%% --------------------------------------
%% --------------------------------------
%% --------------------------------------
%% --------------------------------------
%% --------------------------------------
%% --------------------------------------

end%function
