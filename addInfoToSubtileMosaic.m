function addInfoToSubtileMosaic(subTileDir,dx,outName,varargin)
% addInfoToSubtileMosaic loop through subtiles for info and add to
% existing mosaic file

if exist(outName,'file') ~= 2
    fprintf('File %s not found, returning\n',outName)
    return
end

variableInfo = who('-file', outName);
if ismember('stripList', variableInfo)
    fprintf('File %s already contains stripList var, returning\n',outName)
    return
end

n = find(strcmpi('extent',varargin));
if ~isempty(n)
    x0 = varargin{n+1}(1);
    x1 = varargin{n+1}(2);
    y0 = varargin{n+1}(3);
    y1 = varargin{n+1}(4);
end

fprintf('Indexing subtiles\n')

% make a cellstr of resolved subtile filenames
subTileFiles=dir([subTileDir,'/*_',num2str(dx),'m.mat']);
if isempty(subTileFiles)
    error('No files found matching %s',...
        [subTileDir,'/*_',num2str(dx),'m.mat']);
end

subTileFiles=cellfun( @(x) [subTileDir,'/',x],{subTileFiles.name},...
    'uniformoutput',0);

% Get rows and cols of subtile from file names - assumes the subtile
% names is *_{row}_{col}_{res}.mat
[~,subTileName] = cellfun(@fileparts,subTileFiles,'uniformoutput',0);
subTileName=cellfun(@(x) strsplit(x,'_'),subTileName,'uniformoutput',0);
subTileCol = cellfun(@(x) str2num(x{end-1}),subTileName);
subTileRow = cellfun(@(x) str2num(x{end-2}),subTileName);

Ncols = max(subTileCol);
Nrows = max(subTileRow);

% Get tile projection information, esp. from UTM tile name
%[tileProjName,projection] = getProjName(subTileName{1}{1},projection);

% column-wise subtile number
% 'subTileNum' number DOES NOT MATCH old subtile naming scheme numbers,
% and is for purposes of sorting subtile files only!
subTileNum = sub2ind([Nrows,Ncols],subTileRow,subTileCol);

% % sort subtilefiles by ascending subtile number order
[subTileNum,n] = sort(subTileNum);
subTileCol = subTileCol(n);
subTileRow= subTileRow(n);
subTileFiles = subTileFiles(n);

NsubTileFiles = length(subTileFiles);

% if extent of mosaic not specidied, get from tiles at edges
if ~exist('x0','var')

    % left
    [~,n] = min(subTileCol);
    m = matfile(subTileFiles{n(1)});
    x0 = min(m.x);

    % right
    [~,n] = max(subTileCol);
    m = matfile(subTileFiles{n(1)});
    x1 = max(m.x);

     % bottom
    [~,n] = min(subTileRow);
    m = matfile(subTileFiles{n(1)});
    y0 = min(m.y);

    % top
    [~,n] = max(subTileRow);
    m = matfile(subTileFiles{n(1)});
    y1 = max(m.y);

end

% make a polyshape out of boundary for checking subtile overlap
tilePoly = polyshape([x0 x0 x1 x1]',[y0 y1 y1 y0]');

% build tile coordinate vectors
x = x0:dx:x1;
y = y1:-dx:y0;
  
stripList=[];
parfor filen=1:NsubTileFiles

    fprintf('extracting info from subtile %d of %d: %s\n',filen,NsubTileFiles,subTileFiles{filen})
    
    % get list of variables within this subtile mat file
    mvars  = who('-file',subTileFiles{filen});
    
    if ~any(ismember(mvars,'za_med'))
        fprintf('no za_med found, skipping\n')
        continue
    end

    % Check if subtile overlaps mosaic boundary:
    % open matfile to load coordinate vectors withoutloading arrays
    m=matfile(subTileFiles{filen});

    % make a polyshape out of this subtile boundary
    filenPoly = polyshape([min(m.x) min(m.x) max(m.x) max(m.x)]',...
        [min(m.y) max(m.y) max(m.y) min(m.y)]');

    % test overlap
    if ~overlaps(tilePoly,filenPoly)
        fprintf('subtile out of bounds, skipping\n')
        continue
    end

    % load subtile into structure
    zsub=load(subTileFiles{filen},'x','y');

    % get corner indexes for this grid in the master
    col0  = find(zsub.x(1) == x);
    col1  = find(zsub.x(end) == x);
    row0  = find(zsub.y(1) == y);
    row1  = find(zsub.y(end) == y);

    if isempty(col0) || isempty(col1) || isempty(row0) || isempty(row1)
        fprintf('subtile and tile coordinates dont match, subtile may span boundary,skipping\n')
        continue
    end
    
    % load subtile into structure - iif 2m, need to switch to 10m filename
    % since no filename list in 2m version.
    if any(ismember(mvars,'fileNames'))
        
        zsub=load(strrep(subTileFiles{filen},'_2m.mat','_10m.mat'),'fileNames','dZ');
        zsub.fileNames = zsub.fileNames(~isnan(zsub.dZ));
    
        if isempty(zsub.fileNames)
                fprintf('no filenames returned, skipping\n')
                continue
        end
    
        [~,stripid] =  cellfun(@fileparts,zsub.fileNames,'uniformoutput',0);
        if dx == 2
            stripid= strrep(stripid,'_10m','');
        end
        % stripid =  cellfun(@(x) x(1:47),stripid,'uniformoutput',0);
       
    elseif any(ismember(mvars,'stripIDs'))
        
        zsub=load(subTileFiles{filen},'stripIDs');

        if isempty(zsub.stripIDs)
            fprintf('no stripIDs returned, skipping\n')
            continue
        end
       
        stripid=zsub.stripIDs;
   else
        error('neither stripIDs or fileNames vars in %s',subTileFiles{filen})
    end
       
    stripList = [stripList, stripid];

end

stripList=unique(stripList);
save(outName,'-append','stripList')






