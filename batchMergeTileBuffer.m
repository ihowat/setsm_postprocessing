function batchMergeTileBuffer(fileNames)
% batchMergeTileBuffer: mergeTileBuffer to all neighboring files
%
% batchMergeTileBuffer(fileNames) where fileNames is a cellstr of files to be merged.
fileNames=fileNames(:);

[ntop,nright] = findNeighborTiles(fileNames);

n0 = [nright(:,1);ntop(:,1)];
n1 = [nright(:,2);ntop(:,2)];
    

% run mergeTileBuffer on each pair in list
for i=1:length(n0)
    mergeTileBuffer(fileNames{n0(i)},fileNames{n1(i)},'zOnly'); 
end
    
    
    




