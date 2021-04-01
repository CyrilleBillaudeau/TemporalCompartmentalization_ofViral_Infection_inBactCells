function [tab_ThldInf]=estimateActivationThreshold(lstFolder,nFolder)
%% return in the array with folderTD, and the fluo intensity (min and max) 
%% of the kmeans segmentation (nK=3) of bgd, cell cytoplasm (outside or wo/ SPP1 DNA),
%% SPP1 DNA

tab_ThldInf=zeros(nFolder,7);
for iFolder=1:nFolder
    % load images: binary mask of cells and fluorescent channel used to
    % visualize SPP1 DNA (RFP here)
    disp('===============================================================')
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
        imgRFP=double(imread(imgFilename,2)); % SPP1 DNA
        propROI=regionprops(mskCell,'Area','Centroid','Orientation','BoundingBox','PixelList');
        numROI=numel(propROI);
        figure(16);clf;imagesc(mskCell);cmap=[0,0,0;lines(numROI)];colormap(cmap);
 
        disp('Infected cells identification ...');
        imgInfection=medfilt2(imgRFP,[3 3]);
        
        valBgd=min(imgInfection(mskCell>0));
        imgInfection(mskCell==0)=valBgd; % only for display improvements
        figure(17);clf;imagesc(imgInfection);colormap(hot);
        
        imgInfectionGlobal=kmeans2D(imgInfection,3);close(50);
        tab_ThldInf(iFolder,:)=[iFolder,...
            min(imgInfection(imgInfectionGlobal==1)),max(imgInfection(imgInfectionGlobal==1)),...
            min(imgInfection(imgInfectionGlobal==2)),max(imgInfection(imgInfectionGlobal==2)),...
            min(imgInfection(imgInfectionGlobal==3)),max(imgInfection(imgInfectionGlobal==3))];
        
    end%if
end%for

end%function