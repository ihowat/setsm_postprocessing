function tileNeighborIndexFile=tileNeighborIndex(tileDir)

fileNames = selectTilesForMerging(tileDir);

% find neighbors
[ntop,nright,ntopRight,nbottomRight] = findNeighborTiles(fileNames);

nN = nan(length(fileNames),8);
nN(ntop(:,1),1) = ntop(:,2); %top
nN(ntop(:,2),2) = ntop(:,1); % bottom
nN(nright(:,2),3) = nright(:,1); %left
nN(nright(:,1),4) = nright(:,2); %right
nN(nbottomRight(:,2),5) = nbottomRight(:,1); %top-left
nN(ntopRight(:,1),6) = ntopRight(:,2); %top-right
nN(ntopRight(:,2),7) = ntopRight(:,1); %bottom-left
nN(nbottomRight(:,1),8) = nbottomRight(:,2); %bottom-right

tileNeighborIndexFile=[tileDir,'/tileNeighborIndex.mat'];
save(tileNeighborIndexFile,'fileNames','nN')