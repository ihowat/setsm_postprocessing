function tileNeighborIndexFile=tileNeighborIndex(tileDir,varargin)

n = find(strcmpi('resolution',varargin));
if ~isempty(n)
    resolution = varargin{n+1};
    use_res_in_outname = true;
else
    resolution = '10m';
    use_res_in_outname = false;
end

fileNames = selectTilesForMerging(tileDir,'resolution',resolution);

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

if use_res_in_outname
    tileNeighborIndexFile=[tileDir,'/tileNeighborIndex_',resolution,'.mat'];
else
    tileNeighborIndexFile=[tileDir,'/tileNeighborIndex.mat'];
end
fprintf('writing %s\n', tileNeighborIndexFile)
save(tileNeighborIndexFile,'fileNames','nN')
