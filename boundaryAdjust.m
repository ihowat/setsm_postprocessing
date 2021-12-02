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

%% dz0 calc and write loop - this could be run simultaneously/parallel
i=1;
for i =1:length(fileNames)
    fprintf('Calculating offset grid for %s .... ',fileNames{i})
    
    m0=matfile(fileNames{i});
    if any(strcmp(fields(m0),'dz0'))
            fprintf('dz0 already exists, skipping\n')
            continue
    end
    
    calcdz0(fileNames{i},resizeFraction)
    fprintf('done\n')
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
    
   m0.Properties.Writable = true;
   m0.z = m0.z - m0.dz0;
   m0.adjusted = true;
   fprintf('done\n')
end

function calcdz0(fileName,resizeFraction)

m0=matfile(fileName);
s=whos(m0);
sz=  s(strcmp({s.name},'z')).size;
dz0 = nan(ceil(sz.*resizeFraction));

[path,name,ext]=fileparts(fileName);

if isempty(path)
    path='.';
end

a = strsplit(name,'_');
stringLength=num2str(length(a{1}));
tileRow=str2num(a{1});
tileCol=str2num(a{2});

formatString = ['%s/%0',stringLength,'d_%0',stringLength,'d'];

% above
file1 = [sprintf(formatString,path,tileRow+1,tileCol),...
sprintf('_%s',a{3:end}),ext];
if exist(file1,'file')
    m1=matfile(file1);
    dzbuff = getdzbuff(m0,m1);
    dz0(1,:) = imresize(nanmean(dzbuff),resizeFraction);
else
    dz0(1,:) = 0;
end

% below
file1 = [sprintf(formatString,path,tileRow-1,tileCol),...
sprintf('_%s',a{3:end}),ext];
if exist(file1,'file')
    m1=matfile(file1);
    dzbuff = getdzbuff(m0,m1);
    dz0(end,:) =imresize(nanmean(dzbuff),resizeFraction);
else
    dz0(end,:) = 0;
end

% left
file1 = [sprintf(formatString,path,tileRow,tileCol-1),...
sprintf('_%s',a{3:end}),ext];
if exist(file1,'file')
    m1=matfile(file1);
    dzbuff = getdzbuff(m0,m1);
    dz0(:,1) = imresize(nanmean(dzbuff,2),resizeFraction);
else
    dz0(:,1) = 0;
end

% right
file1 = [sprintf(formatString,path,tileRow,tileCol+1),...
sprintf('_%s',a{3:end}),ext];
if exist(file1,'file')
    m1=matfile(file1);
    dzbuff = getdzbuff(m0,m1);
    dz0(:,end) = imresize(nanmean(dzbuff,2),resizeFraction);
else
    dz0(:,end) = 0;
end

dz0=single(inpaint_nans(double(dz0),2));

dz0 = imresize(dz0,sz);

dz0 = dz0./2;

m0.Properties.Writable = true;
m0.dz0 = dz0;
m0.adjusted = false;

function dzbuff=getdzbuff(m0,m1)

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

        
        