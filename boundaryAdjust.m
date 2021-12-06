function boundaryAdjust(fileNames,varargin)
% boundaryAdjust adjust tile based on offsets with neighbors over buffers
%
% boundaryAdjust(fileNames) where fileNames is a cell of full filenames of
% tile files. Files must be named row_col_*.mat, where the row and col are
% the tile row and column numbers in the tile arrangment (so the tile above
% would be (row+1)_col_*.mat . The tiles above, below, left and right are
% loading and the offset over the buffer is caculated. The values are then
% upscaled by the resizeFraction, which should be ~ the fractional width 
% of the buffer, and added along the edges of  a correction field (dz0) 
% that is then filled  through interpolation. dz0 is generated for all
% tiles in the list and then, afterward, is applied to all tiles in the
% list. A flag, adusted, is added with value false if dz0 has not been
% applied and true if it has.
%
% boundaryAdjust(fileNames,'noApply')only calculates/saves the dz0 field
% but does not apply it to z.

resizeFraction=0.1;

if ~iscell(fileNames)
    fileNames = {fileNames};
end

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

%% dz0 calc and write loop - this could be run simultaneously/parallel
i=1;
for i =1:length(fileNames)
    fprintf('Calculating offset grid for %s .... ',fileNames{i})
    
    m0=matfile(fileNames{i});
    if any(strcmp(fields(m0),'dz0'))
            fprintf('dz0 already exists, skipping\n')
            continue
    end
    
    neighborFiles = cell(8,1);
    neighborFiles(~isnan(nN(i,:))) = fileNames(nN(i,~isnan(nN(i,:))));
    
    calcdz0(fileNames{i},neighborFiles,resizeFraction)
    fprintf('done\n')
    clear  neighborFiles
end

if ~isempty(varargin)
    if any(strcmpi(varargin,'noapply'))
        return
    end
end


%% dz0 apply loop - this could be run simultaneously/parallel
i=1;
for i =1:length(fileNames)
    fprintf('Applying offset to %s ....',fileNames{i})
    
     m0=matfile(fileNames{i});
     
    if any(strcmp(fields(m0),'adjusted'))
        if m0.adjusted 
            fprintf('adjusted flag already true, skipping\n')
            continue
        end
    end
    
    if ~any(strcmp(fields(m0),'dz0'))
            fprintf('dz0 doesnt exist, skipping\n')
            continue
    end
    
    
    s=whos(m0);
    sz=  s(strcmp({s.name},'z')).size;
    
    m0.Properties.Writable = true;
    dz0 = imresize(m0.dz0,sz);
    m0.z = m0.z - dz0;
    m0.adjusted = true;
    clear dz0
    fprintf('done\n')
end

function calcdz0(fileName,neighborFiles,resizeFraction)
% calcdz0 claculate neighbor offset field
%
% calcdz0(fileName,neighborFiles,resizeFraction) calculates the z offset
% over the neighboring files, listed in neighborFiles, and interpolates an
% offset field over z. The neighborFiles list is in order: top, bottom,
% left, right, top-left, top-right, bottom-left, bottom-right, with empty
% cells indicating no file.

m0=matfile(fileName);
s=whos(m0);
sz=  s(strcmp({s.name},'z')).size;
dz0 = nan(ceil(sz.*resizeFraction));

% old name-based neighbor selection - wont work for 2m
% [path,name,ext]=fileparts(fileName);
% 
% if isempty(path)
%     path='.';
% end
% 
% a = strsplit(name,'_');
% stringLength=num2str(length(a{1}));
% tileRow=str2num(a{1});
% tileCol=str2num(a{2});
% 
% formatString = ['%s/%0',stringLength,'d_%0',stringLength,'d'];

% above
% file1 = [sprintf(formatString,path,tileRow+1,tileCol),...
% sprintf('_%s',a{3:end}),ext];
file1 = neighborFiles{1};
if exist(file1,'file')
    m1=matfile(file1);
    dzbuff = getdzbuff(m0,m1);
    dzbuff = nanmean(dzbuff);
    dzbuff = fillmissing(dzbuff,'nearest');
    dz0(1,:) = imresize(dzbuff,resizeFraction);
end

% below
% file1 = [sprintf(formatString,path,tileRow-1,tileCol),...
% sprintf('_%s',a{3:end}),ext];
file1 = neighborFiles{2};
if exist(file1,'file')
    m1=matfile(file1);
    dzbuff = getdzbuff(m0,m1);
    dzbuff = nanmean(dzbuff);
    dzbuff = fillmissing(dzbuff,'nearest');
    dz0(end,:)= imresize(dzbuff,resizeFraction);
end

N0 = uint8(~isnan(dz0));

% left
% file1 = [sprintf(formatString,path,tileRow,tileCol-1),...
% sprintf('_%s',a{3:end}),ext];
file1 = neighborFiles{3};
if exist(file1,'file')
    m1=matfile(file1);
    dzbuff = getdzbuff(m0,m1);
    dzbuff = nanmean(dzbuff,2);
    dzbuff = fillmissing(dzbuff,'nearest');
    dz1= imresize(dzbuff,resizeFraction);
    dz0(:,1)=nanmean([dz1 dz0(:,1)],2);
    N0(:,1) = N0(:,1) + uint8(~isnan(dz1)); 
end

% right
% file1 = [sprintf(formatString,path,tileRow,tileCol+1),...
% sprintf('_%s',a{3:end}),ext];
file1 = neighborFiles{4};
if exist(file1,'file')
    m1=matfile(file1);
    dzbuff = getdzbuff(m0,m1);
    dzbuff = nanmean(dzbuff,2);
    dzbuff = fillmissing(dzbuff,'nearest');
    dz1= imresize(dzbuff,resizeFraction);
    dz0(:,end) = nanmean([dz1 dz0(:,end)],2);
    N0(:,end) = N0(:,end) + uint8(~isnan(dz1)); 
end

% upper-left
% file1 = [sprintf(formatString,path,tileRow+1,tileCol-1),...
% sprintf('_%s',a{3:end}),ext];
file1 = neighborFiles{5};
if exist(file1,'file')
    m1=matfile(file1);
    dzbuff = getdzbuff(m0,m1);
    dzbuff = fillmissing(dzbuff,'nearest');
    dzbuff = fillmissing(dzbuff,'nearest',2);
    dzbuff = imresize(dzbuff,resizeFraction);
    
    szb = size(dzbuff);
    
    dzbuff(2:end,2:end) = NaN;
    
    dz0sub = dz0(1:szb(1),1:szb(2));
    N0sub =  single(N0(1:szb(1),1:szb(2)));
    
    dz0(1:szb(1),1:szb(2)) = ...
        ( (dz0sub.*N0sub) + dzbuff ) ./ (N0sub + 1) ;
end

% upper-right
% file1 = [sprintf(formatString,path,tileRow+1,tileCol+1),...
% sprintf('_%s',a{3:end}),ext];
file1 = neighborFiles{6};
if exist(file1,'file')
    m1=matfile(file1);
    dzbuff = getdzbuff(m0,m1);
    dzbuff = fillmissing(dzbuff,'nearest');
    dzbuff = fillmissing(dzbuff,'nearest',2);
    dzbuff = imresize(dzbuff,resizeFraction);
    
    szb = size(dzbuff);
    
    dzbuff(2:end,1:end-1) = NaN;
    
    dz0sub = dz0(1:szb(1),end-szb(2)+1:end);
    N0sub =  single(N0(1:szb(1),end-szb(2)+1:end));
    
    dz0(1:szb(1),end-szb(2)+1:end) = ...
        ( (dz0sub.*N0sub) + dzbuff ) ./ (N0sub + 1);

end

% lower-left
% file1 = [sprintf(formatString,path,tileRow-1,tileCol-1),...
% sprintf('_%s',a{3:end}),ext];
file1 = neighborFiles{7};
if exist(file1,'file')
    m1=matfile(file1);
    dzbuff = getdzbuff(m0,m1);
    dzbuff = fillmissing(dzbuff,'nearest');
    dzbuff = fillmissing(dzbuff,'nearest',2);
    dzbuff = imresize(dzbuff,resizeFraction);
    
    szb = size(dzbuff);
    
    dzbuff(1:end-1,2:end) = NaN;
    
    dz0sub = dz0(end-szb(1)+1:end,1:szb(2));
    N0sub =  single(N0(end-szb(1)+1:end,1:szb(2)));
    
    dz0(end-szb(1)+1:end,1:szb(2)) = ...
        ( (dz0sub.*N0sub) + dzbuff ) ./ (N0sub + 1) ;
end

% lower-right
% file1 = [sprintf(formatString,path,tileRow-1,tileCol+1),...
% sprintf('_%s',a{3:end}),ext];
file1 = neighborFiles{8};
if exist(file1,'file')
    m1=matfile(file1);
    dzbuff = getdzbuff(m0,m1);
    dzbuff = fillmissing(dzbuff,'nearest');
    dzbuff = fillmissing(dzbuff,'nearest',2);
    dzbuff = imresize(dzbuff,resizeFraction);
    
    szb = size(dzbuff);
    
    dzbuff(1:end-1,1:end-1) = NaN;
    
    dz0sub = dz0(end-szb(1)+1:end,end-szb(2)+1:end);
    N0sub =  single(N0(end-szb(1)+1:end,end-szb(2)+1:end));
    
    dz0(end-szb(1)+1:end,end-szb(2)+1:end) = ...
        ( (dz0sub.*N0sub) + dzbuff ) ./ (N0sub + 1);
end

dz0=single(inpaint_nans(double(dz0),2));

dz0 = dz0./2;

m0.Properties.Writable = true;
m0.dz0 = dz0;
m0.adjusted = false;

function dzbuff=getdzbuff(m0,m1)
% calulate the offset over the tile buffer

c0 = m0.x >= min(m1.x) & m0.x <= max(m1.x);
r0 = m0.y >= min(m1.y) & m0.y <= max(m1.y);
c0 = [find(c0,1,'first'),find(c0,1,'last')];
r0 = [find(r0,1,'first'),find(r0,1,'last')];

if isempty(c0) || isempty(r0); error('no overlap between tiles'); end

c1 = m1.x >= min(m0.x) & m1.x <= max(m0.x);
r1 = m1.y >= min(m0.y) & m1.y <= max(m0.y);
c1 = [find(c1,1,'first'),find(c1,1,'last')];
r1 = [find(r1,1,'first'),find(r1,1,'last')];

dzbuff = m0.z(r0(1):r0(2),c0(1):c0(2)) - m1.z(r1(1):r1(2),c1(1):c1(2));

s=whos(m0);
if any(strcmp({s.name},'land'))
    land = m0.land(r0(1):r0(2),c0(1):c0(2));
    dzbuff(~land) = NaN;
end

s=whos(m1);
if any(strcmp({s.name},'land'))
    land = m1.land(r0(1):r0(2),c0(1):c0(2));
    dzbuff(~land) = NaN;
end

        
        