function cellSummary=extractTotInt_InfectedCells(imgRFP,mskCell,cellSummary)
%% Measure the average and std intensity in infected bacteria in the fluorescence 
%% channel of SPP DNA (RFP here) (mono- and multi-infection).
%% Results are written in cellSummary array (columns 4 & 5).

disp('Intensity measurement of infected cells ...');
nROI=size(cellSummary,1);
for iROI=1:nROI
    if (cellSummary(iROI,3)>0)
        val=imgRFP(mskCell==iROI);
        cellSummary(iROI,4:5)=[nanmean(val), nanstd(val)];
    else
        cellSummary(iROI,4:5)=nan(1,2);
    end
end%for
disp('Intensity measurement of infected cells done');
end%function