function batchMergeTileBuffer(tileNeighborIndexFile)
% batchMergeTileBuffer: mergeTileBuffer to all neighboring files
%
% batchMergeTileBuffer((tileNeighborIndexFile) where tileNeighborIndexFile 
% contains the list of fileNames to be merged (fileNames) and the tile
% neighbor index array nN=[indTop, indBottom, indLeft, indRight], created
% by tileNeighborIndex 

load(tileNeighborIndexFile,'fileNames','nN');

% find files where indices of top and right files exist
indTopExists = find(~isnan(nN(:,1)));
indRightExists = find(~isnan(nN(:,4)));

% concatenate vectors of indices of bottom/top and left/right pairs
n0 =[indTopExists; indRightExists]; %[bottom file0; left file0] 
n1 = [nN(indTopExists,1); nN(indRightExists,4)]; %[top file1; right file1]

% run mergeTileBuffer on each pair in list
for i=1:length(n0)
    mergeTileBuffer(fileNames{n0(i)},fileNames{n1(i)}); 
end
    
    
    




