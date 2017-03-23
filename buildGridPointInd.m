function [meta] = buildGridPointInd(meta,x,y,N)
% buildGridPointInd build grid point coverage index for each strip
%
% meta = buildGridPointInd(meta,x,y) make a cell for each file containing
% the col-wise indices of the data points on the tile grid. This is used
% to search for new or existing data on the grid within each search
% iteration in strips2tile.

% first make sure does not exist in case we're saving this info
if ~isfield(meta,'gridPointInd')
    
    fprintf('calculating tile pixel coverage for each strip\n');
    
    % initialize output field  
    meta.gridPointInd=cell(size(meta.f));
    
    % can be slow, so we'll use a counter
    count = 0;
    
    % file loop
    for i=1:length(meta.f)
        
        % counter
        if i>1 for p=1:count fprintf('\b'); end; %delete line before
            count = fprintf('strip %d of %d',i,size(meta.f,1));
        end
        
        % locate grid pixels within footprint polygon
        BW = roipoly(x, y, N, meta.x{i}, meta.y{i});
        
        % if mask exists, apply it
        if isfield(meta,'maskPolyx')
            mask=[meta.maskPolyx{i}(:) meta.maskPolyy{i}(:)];
            
            for j=1:size(mask,1)
                if ~isempty(mask{j,1}) && ~isempty(mask{j,2})
                    BW(roipoly(x,y,N,mask{j,1},mask{j,2}))=0;
                end
            end
            
        end
   
        % convert BW mask to col-wise indices and save to cell
        meta.gridPointInd{i}=find(BW);
        
        % get rid of mask
        clear BW
    end
    
    % clear counter line
    for p=1:count fprintf('\b'); end;
    
    % make another field with the number of data points
    meta.gridPointN=cellfun(@length,meta.gridPointInd);

end
