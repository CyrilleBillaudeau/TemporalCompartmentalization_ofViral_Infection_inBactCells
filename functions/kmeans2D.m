function img_Cluster=kmeans2D(img,nK)
%% Segmentation with k-means

% transform 2D matrix into 1D tab for kmeans analysis
nrows = size(img,1);
ncols = size(img,2);
img_1D=img(:);

%nK=2; % number of clusters
doKemans=1;
while doKemans
    try
        [IDX, meanCl] = kmeans(img_1D, nK,'distance','sqEuclidean','Replicates',3,'start','cluster','emptyaction','singleton');
        doKemans=0;
    catch
        disp('warning error in kmeans');
        disp('kmeans is repeated again');
        doKemans=1;
    end
end
%disp(meanCl);

% cluster sorting with increasing mean values of intensity
%ProjZ_Cluster=zeros(nrows,ncols);
[~,isort]=sort(meanCl);%[meanCl,isort,z];
IDX_sort=zeros(nrows*ncols,1);

for iCl=1:nK
    IDX_sort(find(IDX==isort(iCl)))=iCl;
end

% transform the results of kmeans (1D tab) into a 2D matrix used to get the
% binary mask
IDX=IDX_sort;
img_Cluster = reshape(IDX,nrows,ncols);

figure(50);%%set(70,'WindowStyle', 'docked');
image((1:ncols),(1:nrows),img_Cluster);colormap(gray(nK));% spatial distribution of the clusters
axis equal;colormap(gray);
end%function
