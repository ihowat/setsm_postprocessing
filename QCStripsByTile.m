function QCStripsByTile(regionNum,varargin)

%% Argins
%regionNum='02'; % region number
tilefile  = 'PGC_Imagery_Mosaic_Tiles_Antarctica.mat'; %PGC/NGA Tile definition file. 
arcdemfile= 'rema_tiles.mat'; % lists which tiles go to which regions
%dbasefile = 'rema_strips_8m_wqc_cs2bias.mat';
dbasefile = 'rema_strips_2m_wbias.mat';
startfrom = 1;
minN = 500;
minArea = 500;

for i=1:2:length(varargin)
    eval([varargin{i},'=',num2str(varargin{i+1}),';']);
end

tileDir= dir(['/data4/REMA/region_',regionNum,'*']);
%tileDir= ['/data4/REMA/',tileDir(1).name,'/mosaic_reg_qc_feather2/40m/'];
tileDir= ['/data4/REMA/',tileDir(1).name,'/mosaic2m/40m/'];


%Get Arctic Tile Defs
tiles=load(tilefile);
a=load(arcdemfile);

%remove duplicated tile entries keeping just first occurence
[~,n]  = unique(a.tileName,'stable');
a = structfun(@(x) ( x(n) ), a, 'UniformOutput', false);

% get index of this region number
n=a.regionNum==str2num(regionNum);

% check to make sure we find sum
if ~(any(n)); error('no tiles matched for this region number'); end

% match region number to overlapping tiles
[~,n]=intersect(tiles.I,a.tileName(n));

% crop tile structure to overlapping tiles
tiles = structfun(@(x) ( x(n) ), tiles, 'UniformOutput', false);

% load database structure
meta=load(dbasefile);

%check meta file for required fields
flds = fields(meta);

if ~any(strcmp(flds,'avg_rmse')); error('meta stucture missing avg_rmse field \n'); end
if ~any(strcmp(flds,'xmax')); error('meta stucture missing xmax field \n'); end
if ~any(strcmp(flds,'ymax')); error('meta stucture missing ymax field \n'); end
if ~any(strcmp(flds,'xmin')); error('meta stucture missing ymin field \n'); end
if ~any(strcmp(flds,'ymin')); error('meta stucture missing ymin field \n'); end
if ~any(strcmp(flds,'x')); error('meta stucture missing x field \n'); end
if ~any(strcmp(flds,'y')); error('meta stucture missing y field \n'); end
if ~any(strcmp(flds,'f')); error('meta stucture missing f field \n'); end
if ~any(strcmp(flds,'f')); error('meta stucture missing f field \n'); end
if ~any(strcmp(flds,'sigma_all')); error('meta stucture missing sigma_all field \n'); end
if ~any(strcmp(flds,'sigma_1yr')); error('meta stucture missing sigma_1yr field \n'); end


% select the whichever registration has the better sigma_bias (all or 1 yr)
meta.sigma = nanmin([meta.sigma_all(:)';meta.sigma_1yr(:)'])';
meta.sigma(meta.sigma > 1) = NaN;

meta.avg_rmse(meta.avg_rmse == 0) = NaN;

for i=startfrom:length(tiles.I)
    fprintf(' \n')
    fprintf('Working tile %d of %d: %s \n',i,length(tiles.I),tiles.I{i});
    tile = structfun(@(x) ( x(i) ), tiles, 'UniformOutput', false);
    
    if exist([tileDir,tile.I{1},'_40m_dem.mat'],'file')
        load([tileDir,tile.I{1},'_40m_dem.mat'],'N');
        
        if sum(N(:))./numel(N) > 0.999
            
            fprintf('tile coverage complete, skipping\n')
            
            clear N
            
            continue
            
        end
        
    end
    
    
    qctile(tile,meta,minN,minArea);
    
end

function qctile(tiles,meta,minN,minArea)

%% Spatial coverage search

% make tile boundary polygon

tilevx = [tiles.x0;tiles.x0;tiles.x1;tiles.x1;tiles.x0];
tilevy = [tiles.y0;tiles.y1;tiles.y1;tiles.y0;tiles.y0];

% quick search: find strips within range of this tile. This does not
% account for background area of around strips but just pairs them down to
% speed the poly intersection loop

n = meta.xmax > tiles.x0 &  meta.xmin < tiles.x1 & ...
    meta.ymax > tiles.y0 &  meta.ymin < tiles.y1;

if ~any(n); fprintf('no strip overlap\n'); return; end

% get polys
meta = structfun(@(x) ( x(n,:) ), meta, 'UniformOutput', false);

% search for all strips overlapping this tile
in=zeros(size(meta.f));

for i=1:length(in)
    in(i) = any(inpolygon(meta.x{i},meta.y{i},tilevx,tilevy)) | ...
        any(inpolygon(tilevx,tilevy,meta.x{i},meta.y{i}));
end

if ~any(in); fprintf('no strip overlap\n'); return; end

% crop meta data struct to only overlapping strips
meta = structfun(@(x) ( x(logical(in),:) ), meta, 'UniformOutput', false);

fprintf('%d files overlapping this tile, ',sum(in));

% add existing qc data
meta = addQC2Meta(meta);

fprintf('%d files with existing qc\n',sum(meta.qc ~= 0));


% remove already added entries  from metadata
if any(meta.qc == 5)
    fprintf('removing %d files with qc flag=5\n',sum(meta.qc == 5));
    meta = structfun(@(x) ( x(meta.qc ~= 5 ,:) ), meta, 'UniformOutput', false);
    if isempty(meta.f); fprintf('all strips removed, returning \n'); return; end
end

% build grid
res = 40;
buff = 0;

x = tiles.x0-buff*res: res:tiles.x1+buff*res;
y = tiles.y1+buff*res:-res:tiles.y0-buff*res;
y = y(:);

N= zeros(length(y),length(x),'uint8');


if isfield(tiles,'coastline')
    fprintf('applying coastline, ');
    
    
    A = false(size(N));
    i=1;
    for i=1:length(tiles.coastline{1})
        if  isempty(tiles.coastline{1}{i}); continue; end
        A(roipoly(x,y,N,tiles.coastline{1}{i}(1,:),...
            tiles.coastline{1}{i}(2,:))) = true;
    end
    
    
    percent_filled = 100*sum(~A(:))./numel(A);
    
    if percent_filled >  0.2  
        N(~A) = 1; 
    else
        percent_filled =0;
    end
    clear A
    
     fprintf('%.2f%% filled as water\n',percent_filled);

    if percent_filled == 100; fprintf('returning \n'); return; end
end

%% Build Grid Point Index Field
% make a cell for each file containing the col-wise indices of the data
% points on the tile grid. This is used search for new or existing data on
% the grid within each search iteration.

fprintf('calculating tile pixel coverage for each strip, ');

% initialize output field
meta.gridPointInd=cell(size(meta.f));

% can be slow, so we'll use a counter
count = 0;

% already qc'd files to add to N
add2N = meta.qc > 0 & meta.qc < 4;

fprintf('%d already qc-passed files will be added \n',sum(add2N))

if ~any(~add2N); fprintf('all files already passed qc, returning \n'); return; end

% file loop
for i=1:length(meta.f)
    
    % counter
    if i>1 for p=1:count fprintf('\b'); end; %delete line before
        count = fprintf('strip %d of %d',i,size(meta.f,1));
    end
    
    % locate grid pixels within footprint polygon
    BW = roipoly(x, y, N, meta.x{i}, meta.y{i});
    
    %if mask data exists, apply it
    if meta.qc(i) == 3
        j=1;
        for j=1:length(meta.maskPolyx{i})
            BW(roipoly(x,y,BW,...
                meta.maskPolyx{i}{j},meta.maskPolyy{i}{j}))=0;
        end
    end
    
    % add if already qc'd
    if add2N(i); N(BW) = 1; continue; end
    
    % convert BW mask to col-wise indices and save to cell
    meta.gridPointInd{i}=find(BW);
    
    % get rid of mask
    clear BW
end

% clear counter line
for p=1:count fprintf('\b'); end;

% check if filled
percent_filled = 100*sum(N(:))./numel(N);
fprintf('%.2f%% filled\n',percent_filled);
if percent_filled == 100; fprintf('returning \n'); return; end

% remove already added entries  from metadata
meta = structfun(@(x) ( x(~add2N ,:) ), meta, 'UniformOutput', false);

% make another field with the number of data points
meta.gridPointN=cellfun(@length,meta.gridPointInd);


%remove strips below a minimum size
stripArea = nan(size(meta.f));
for i=1:length(meta.f)
    stripArea(i) = polyarea(meta.x{i}, meta.y{i})./1000^2;
end

fprintf('removing %d strips smaller than %.1f km^2\n',sum(stripArea < minArea),minArea);
meta = structfun(@(x) ( x(stripArea >= minArea ,:) ), meta, 'UniformOutput', false);



%% coverage test loop
if ~isfield(meta,'rmse'); meta.rmse=nan(size(meta.f));end
if ~isfield(meta,'dtrans'); meta.dtrans=zeros(length(meta.f),3); end
if ~isfield(meta,'overlap'); meta.overlap=zeros(size(meta.f)); end % number existing data points overlapping file
%%
recountFlag=true;
skipn = 0;
while length(meta.f) >= 1
    
    percent_filled = 100*sum(N(:))./numel(N);
    
    if percent_filled == 100; fprintf('100%% filled returning \n',percent_filled); return; end
    
    fprintf('%.2f%% filled\n', percent_filled);
    
    % loop through files in boundary and find # of new/existing points
    if recountFlag
        for i=1:length(meta.f)
            
            % subset of N at data points in this file
            Nsub=N(meta.gridPointInd{i});
            
            % count existing data points
            meta.overlap(i)=sum(Nsub);
            
            % count new data points
            meta.gridPointN(i) = sum(~Nsub);
            
        end
    end
    
    recountFlag=false;
    
    % remove rendundant files (with already 100% coverage)
    redundantFlag = meta.gridPointN < minN;
    
    if any(redundantFlag)

        % remove redundant files from lists
        meta = structfun(@(x) ( x(~redundantFlag,:) ), meta, 'UniformOutput', false);
        
        fprintf('%d redundant files (N < %d) removed\n',sum(redundantFlag),minN);
        
        if isempty(meta.f); fprintf('all strips removed, returning \n'); return; end
    end

    A = nansum([100.*meta.gridPointN./numel(N),1./(meta.avg_rmse.^2)],2);
 
    [~,n]=sort(A,'descend');
    
    
    if skipn < 0; skipn = length(n)-1; end
    
    % skip if skipped on last iteration
    if length(n) >= 1+skipn
        n = n(1+skipn);
    else
        n=n(1);
        skipn = 0;
    end
    
    fprintf('%d of %d strips remaining\n',skipn+1,length(meta.f));

    fileName=strrep(meta.f{n},'meta.txt','dem_browse.tif');
    
    
    fprintf('%s\n',fileName);
    fprintf('%d new pointsm, gcp sigma=%.2f, mean coreg RMSE=%.2f, max coreg RMSE=%.2f \n',...
        meta.gridPointN(n),meta.sigma(n),meta.avg_rmse(n), meta.max_rmse(n));
    
    % localfileName = strrep(fileName,'/data4/REMA','~/rema8mStripBrowse');
    localfileName = fileName;
    
    I=readGeotiff(localfileName);
     
    Ni=interp2(x,y(:),N,I.x,I.y(:),'*nearest');
    
    imagesc(I.x,I.y,I.z,'alphadata',single(I.z ~= 0))
    set(gca,'color','r')
    axis xy  equal tight;
    colormap gray;
    hold on;
    imagesc(I.x,I.y,Ni,'alphadata',single(Ni).*.5)
    

    if isfield(tiles,'coastline')
        i=1;
        for i=1:length(tiles.coastline{1})
            if  isempty(tiles.coastline{1}{i}); continue; end
            plot(tiles.coastline{1}{i}(1,:),...
                tiles.coastline{1}{i}(2,:),'b','linewidth',1)
        end
    end
    
    plot([tiles.x0,tiles.x0,tiles.x1,tiles.x1,tiles.x0], [tiles.y0,tiles.y1,tiles.y1,tiles.y0,tiles.y0],'w','linewidth',2)
   
    set(gca,'xlim',[min(I.x)-500 max(I.x)+500],'ylim',[min(I.y)-500 max(I.y)+500]);
    
    
    %set(gcf,'units','normalized');
    %set(gcf,'position',[0.01,0.01,.35,.9])
    
    
    qc=load([fileparts(fileName),'/qc.mat']);
    
    % qc.fileNames=strrep(qc.fileNames,'/data2','/data3');
    
    [~,IA]=intersect(qc.fileNames, fileName);
    
    if isempty(IA) 
        error('this file name not matched in the qc.mat, probably need to upadate it.'); 
    end
    
    if qc.flag(IA) ~= 4 && qc.flag(IA) ~= 0
        
        fprintf('flag previoulsy changed to %d, applying\n',qc.flag(IA))
        
    else
        
        j=1;
        while j
            
            try
                
                flag=input('Enter quality flag: 0=skip,9=back, 1=good, 2=partial, 3=manual edit, 4=poor, 5=unuseable, 6=quit\n');
                
                
                if ~isempty(flag)
                    if isnumeric(flag)
                        if flag == 0 || flag == 9 || flag == 1 || flag == 2 || flag == 3 || flag == 4 || flag == 5 || flag == 6
                            break;
                        end
                    end
                end
            catch
                fprintf('%d not recogized, try again\n',flag);
            end
            
            if iscell(flag); flag=flag{1}; end
            
            fprintf('%d not recogized, try again\n',flag);
            
        end
        
        
        if flag == 6; clf; return; end
        
        
        if flag == 0;  skipn = skipn+1;  clf; continue; end
        
        if flag == 9;  skipn = skipn-1;  clf; continue; end
        
        
        qc.flag(IA)=flag;
        
        if flag == 1 || flag == 2
            qc.x{IA}=cell(1); qc.y{IA}=cell(1);
            
        end
        
        if flag == 3
            fprintf('entering manual edit mode\n')
            
            j=1;
            while j
                
                [~,qc.x{IA}{j},qc.y{IA}{j}] =  roipoly;
                
                plot(qc.x{IA}{j},qc.y{IA}{j},'g','linewidth',2)
                
                
                while j
                    s=input('continue editing this image? (y/n)\n','s');
                    
                    if ~strcmpi(s,'y') && ~strcmpi(s,'n')
                        fprintf('%s not recogized, try again\n',s);
                    else
                        break
                    end
                end
                
                if strcmpi(s,'n'); break; end
                
                j=j+1;
            end
        end
        
        save([fileparts(fileName),'/qc.mat'],'-struct','qc');
        
    end
    
    
    clf
    
    if qc.flag(IA) > 0 && qc.flag(IA) < 4
        
        M = I.z ~=0;
        
        j=1;
        for j=1:length(qc.x{IA})
            M(roipoly(I.x,I.y,M,qc.x{IA}{j},qc.y{IA}{j}))=0;
        end
        
        M=interp2(I.x,I.y(:),single(M),x,y(:),'*nearest');
        
        N(M == 1) = 1;
        
        recountFlag=true;
        
        % remove this file from the meta struct
        in=1:length(meta.f); in(n)=[];
        meta = structfun(@(x) ( x(in,:) ), meta, 'UniformOutput', false);
        skipn = 0;
        
    elseif qc.flag(IA) == 4 || qc.flag(IA) == 5
        
        % remove this file from the meta struct
        in=1:length(meta.f); in(n)=[];
        meta = structfun(@(x) ( x(in,:) ), meta, 'UniformOutput', false);
        
        
    end
    
end


%