function boundaryAdjustCalc(fileName,neighborFiles,varargin)
% boundaryAdjustCalc claculate neighbor offset field
%
% boundaryAdjustCalc(fileName,neighborFiles,resizeFraction) calculates the z offset
% over the neighboring files, listed in neighborFiles, and interpolates an
% offset field over z. The neighborFiles list is in order: top, bottom,
% left, right, top-left, top-right, bottom-left, bottom-right, with empty
% cells indicating no file.

m0=matfile(fileName);

resizeFraction = 0.1;
%% standard resizeFraction is for 10m, scale for 2m and other res
%resizeFraction_10m = 0.1;
%res = m0.x(1,2) - m0.x(1,1);
%resizeFraction = min(1.0, resizeFraction_10m * (res/10));

n=find(strcmpi(varargin,'resizeFraction'));
if ~isempty(n)
    resizeFraction=varargin{n+1};
end

if any(strcmpi(fields(m0),'adjusted')) 
    if m0.adjusted == 1
        fprintf('adjustment applied, undo adjustment before calculating new, skipping, ')
        return
    end
end

sz= size(m0,'z');
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
    [dzbuff,xbuff,ybuff]  = getdzbuff(m0,m1);
    dzbuff = nanmean(dzbuff);
    
    if length(xbuff) < sz(2)
        [~,IA,IB]=intersect(m0.x,xbuff);
        dzbuff1=nan(1,sz(2));
        dzbuff1(IA) = dzbuff(IB);
        dzbuff=dzbuff1;
    end
    
    dzbuff = fillmissing(dzbuff,'nearest');
    dz0(1,:) = imresize(dzbuff,resizeFraction);
end

% below
% file1 = [sprintf(formatString,path,tileRow-1,tileCol),...
% sprintf('_%s',a{3:end}),ext];
file1 = neighborFiles{2};
if exist(file1,'file')
    m1=matfile(file1);
    [dzbuff,xbuff,ybuff]  = getdzbuff(m0,m1);
    dzbuff = nanmean(dzbuff);
    
    if length(xbuff) < sz(2)
        [~,IA,IB]=intersect(m0.x,xbuff);
        dzbuff1=nan(1,sz(2));
        dzbuff1(IA) = dzbuff(IB);
        dzbuff=dzbuff1;
    end
    
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
    [dzbuff,xbuff,ybuff] = getdzbuff(m0,m1);
    dzbuff = nanmean(dzbuff,2);
    
    if length(ybuff) < sz(1)
        [~,IA,IB]=intersect(m0.y,ybuff);
        dzbuff1=nan(sz(1),1);
        dzbuff1(IA) = dzbuff(IB);
        dzbuff=dzbuff1;
    end
    
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
    [dzbuff,xbuff,ybuff]  = getdzbuff(m0,m1);
    dzbuff = nanmean(dzbuff,2);
    
    if length(ybuff) < sz(1)
        [~,IA,IB]=intersect(m0.y,ybuff);
        dzbuff1=nan(sz(1),1);
        dzbuff1(IA) = dzbuff(IB);
        dzbuff=dzbuff1;
    end
    
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
    [dzbuff,xbuff,ybuff]  = getdzbuff(m0,m1);
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
    [dzbuff,xbuff,ybuff]  = getdzbuff(m0,m1);
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
    [dzbuff,xbuff,ybuff]  = getdzbuff(m0,m1);
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
    [dzbuff,xbuff,ybuff]  = getdzbuff(m0,m1);
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

function [dzbuff,x,y]=getdzbuff(m0,m1)
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

dz1 = zeros(length(r1(1):r1(2)),length(c1(1):c1(2)));
if any(strcmpi(fields(m1),'adjusted'))
    if m1.adjusted == 1
        dz1=imresize(m1.dz0,size(m1,'z'));
        dz1=dz1(r1(1):r1(2),c1(1):c1(2));
    end
end

dzbuff = m0.z(r0(1):r0(2),c0(1):c0(2)) - (m1.z(r1(1):r1(2),c1(1):c1(2)) + dz1);

if any(strcmp(fields(m0),'land'))
    land = m0.land(r0(1):r0(2),c0(1):c0(2));
    dzbuff(~land) = NaN;
end

if any(strcmp(fields(m1),'land'))
    land = m1.land(r1(1):r1(2),c1(1):c1(2));
    dzbuff(~land) = NaN;
end

x=m0.x(1,c0(1):c0(2));
y=m0.y(1,r0(1):r0(2));
