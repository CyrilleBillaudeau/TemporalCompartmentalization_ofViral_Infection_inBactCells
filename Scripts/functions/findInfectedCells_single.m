function [tabInfectedCells,imgInfectionGlobal,imgInfectionGlobalHigh,propROI_infection]=findInfectedCells_single(imgRFP,mskCell,limAreaInf)
%% find infected cells using kmeans with RFP fluo (SPP1 DNA)
%% Outputs:
%% tabInfectedCells: number of SPP1 DNA infection per cell
%% imgInfectionGlobal: resulting image of kmeans with 3 classes (1=bgd, 2=cell cytoplasm wo/ infection, 3=SPP1 DNA area)
%% imgInfectionGlobalHigh: SPP1 DNA binary mask
%% propROI_infection: SPP1 DNA region properties ('Area','BoundingBox','Centroid');

disp('Infected cells identification ...');
imgInfection=medfilt2(imgRFP,[3 3]);
valBgd=min(imgInfection(mskCell>0));
imgInfection(mskCell==0)=valBgd; % display improvements
figure(17);clf;imagesc(imgInfection);colormap(hot);

imgInfectionGlobal=kmeans2D(imgInfection,3);close(50);

figure(18);clf;imagesc(imgInfectionGlobal);colormap(gray);

% extract region properties of all SPP1 DNA
imgInfectionGlobalHigh=bwlabel(imgInfectionGlobal==3);
propROI_infection=regionprops(imgInfectionGlobalHigh,'Area','BoundingBox','Centroid');
nROI_inf=numel(propROI_infection);

% remove infection area with area below limAreaInf
for iROI_inf=1:nROI_inf
    if (propROI_infection(iROI_inf).Area < limAreaInf)
        imgInfectionGlobal(imgInfectionGlobalHigh==iROI_inf)=2;
        %disp(strcat(['small :',num2str(iROI_inf)]))
    end
end
nROI_inf0=nROI_inf;

% extract region properties of all SPP1 DNA (after removing small ones)
imgInfectionGlobalHigh=bwlabel(imgInfectionGlobal==3);
propROI_infection=regionprops(imgInfectionGlobalHigh,'Area','BoundingBox','Centroid');
nROI_inf=numel(propROI_infection);
disp(strcat(['Found ',num2str(nROI_inf),' SPP1 DNA after rejecting small areas (',num2str(nROI_inf0-nROI_inf),' were refejected)']));

% quantify SPP1 DNA region inside bacteria to get non-infected, mono- and multi- infected cells
nROI=max(mskCell(:));
tabInfectedCells=zeros(nROI,1);
for iROI=1:nROI
    % isolate SPP1 DNA regions inside current cell (ROI)
    imgInfectionGlobal_ROI=imgInfectionGlobal;
    imgInfectionGlobal_ROI(mskCell~=iROI)=0;
    imgInfectionGlobal_ROI(imgInfectionGlobal_ROI<3)=0;

    % count SPP1 DNA area in current cell
    imgInfectionGlobal_ROI=bwlabel(imgInfectionGlobal_ROI);
    %figure(901);clf;imagesc(imgInfectionGlobal_ROI);colormap(gray);
    tabInfectedCells(iROI)=max(imgInfectionGlobal_ROI(:));
end

disp('Infected cells identification done');
end%function