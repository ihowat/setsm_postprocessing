function [] = QCStripsByTile_above(regionNum,varargin)

%% Argins
%regionNum='02'; % region number
tilefile  = 'V:/pgc/data/scratch/claire/repos/setsm_postprocessing_pgc/PGC_Imagery_Mosaic_Tiles_Above_nocoast.mat'; %PGC/NGA Tile definition file, required
arcdemfile= 'V:/pgc/data/scratch/claire/repos/setsm_postprocessing_pgc/ABoVE_tiles.mat'; % lists which tiles go to which regions, required
dbasefile = 'V:/pgc/data/scratch/claire/repos/setsm_postprocessing_pgc/aboveDEMdatabase_2m.mat'; % database file
changePath= 'V:/pgc'; %if set, will change the path to the REMA directory from what's in the database file. set to [] if none.

% if an older set of mosaic files already exist, we can speed things up by
% check to see if they already have 100% coverage - will skip if do. Leave
% empty if none.
tileDir= '/mnt/pgc/data/elev/dem/setsm/aboveDEM/mosaic/2m_v2/';

startfrom = '1';
minN = 500;
minArea = 35;

for i=1:2:length(varargin)
    eval([varargin{i},'=''',(varargin{i+1}),''';']);
end
startfrom=str2num(startfrom);

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
fprintf('Loading db\n');
meta=load(dbasefile);

% check for region field
if ~isfield(meta,'region')
    meta.region=cell(size(meta.f));
    i=1;
    for i=1:length(meta.f); 
            meta.region{i} = fileparts(meta.f{i});
    end
end

% alter paths in database if set
if ~isempty(changePath)
    meta.f = strrep(meta.f,'/mnt/pgc',changePath);
    meta.f = strrep(meta.f,'/','\');    
    meta.region = strrep(meta.region,'/mnt/pgc',changePath);
    meta.region = strrep(meta.region,'/','\');
end

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

% select the whichever registration has the better sigma_bias (all or 1 yr)
if isfield(meta,'sigma_all') &&  isfield(meta,'sigma_1yr') &&  ~isfield(meta,'sigma') 
    meta.sigma = nanmin([meta.sigma_all(:)';meta.sigma_1yr(:)'])';
end

% if no ground control error field, just set to nan
if ~isfield(meta,'sigma')
    meta.sigma  = nan(size(meta.f));
end

% if ground control error > 1, set to NaN
meta.sigma(meta.sigma > 1) = NaN;
meta.avg_rmse(meta.avg_rmse == 0) = NaN;

for i=startfrom:length(tiles.I)
    fprintf('\n')
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
    
    
    qctile(tile,meta,minN,minArea,changePath);
    
end
end

function qctile(tiles,meta,minN,minArea,changePath)

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
meta = addQC2Meta(meta,changePath);

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
    rawfileName=strrep(meta.f{n},'meta.txt','dem_small.tif');
    
    
    fprintf('%s\n',fileName);
    fprintf('%d new pointsm, gcp sigma=%.2f, mean coreg RMSE=%.2f, max coreg RMSE=%.2f \n',...
        meta.gridPointN(n),meta.sigma(n),meta.avg_rmse(n), meta.max_rmse(n));
    
    % localfileName = strrep(fileName,'/data4/REMA','~/rema8mStripBrowse');
    localfileName = fileName;
    
    I=readGeotiff(localfileName);
     
%     Ni=interp2(x,y(:),N,I.x,I.y(:),'*nearest');
    
    
    % Make square image to show in the figure window.
    X_fig = I.x;
    Y_fig = I.y;
    xsize = max(size(I.x));
    ysize = max(size(I.y));
    max_dim = max(xsize, ysize);
    if max_dim == ysize
        dx = X_fig(2) - X_fig(1);
        x_border = (max_dim - xsize) / 2;
        x_lft = [(X_fig(1)-floor(x_border)*dx):dx:(X_fig(1)-dx)];
        x_rgt = [(X_fig(end)+dx):dx:(X_fig(end)+ceil(x_border)*dx)];
        X_fig = [x_lft, X_fig, x_rgt];
    else
        dy = Y_fig(2) - Y_fig(1);
        y_border = (max_dim - ysize) / 2;
        y_lft = [(Y_fig(1)-floor(y_border)*dy):dy:(Y_fig(1)-dy)];
        y_rgt = [(Y_fig(end)+dy):dy:(Y_fig(end)+ceil(y_border)*dy)];
        Y_fig = [y_lft, Y_fig, y_rgt];
    end
    
    Ni=interp2(x,y(:),N,X_fig,Y_fig(:),'*nearest');
    
    Z_fig = zeros(max(size(Y_fig)), max(size(X_fig)));
    Ni_fig = zeros(max(size(Y_fig)), max(size(X_fig)));
    
    Z_c0 = find(I.x(1) == X_fig);
    Z_c1 = find(I.x(end) == X_fig);
    Z_r0 = find(I.y(1) == Y_fig);
    Z_r1 = find(I.y(end) == Y_fig);
    
    Ni_c0 = find(x(1) == X_fig);
    if isempty(Ni_c0)
        Ni_c0 = 1;
    end
    Ni_c1 = find(x(end) == X_fig);
    if isempty(Ni_c1)
        Ni_c1 = length(X_fig);
    end
    Ni_r0 = find(y(1) == Y_fig);
    if isempty(Ni_r0)
        Ni_r0 = 1;
    end
    Ni_r1 = find(y(end) == Y_fig);
    if isempty(Ni_r1)
        Ni_r1 = length(Y_fig);
    end
    
    Z_fig(Z_r0:Z_r1, Z_c0:Z_c1) = I.z;
    Ni_fig(Ni_r0:Ni_r1, Ni_c0:Ni_c1) = Ni(Ni_r0:Ni_r1, Ni_c0:Ni_c1);
    
    
    fig = figure;
    imagesc(X_fig,Y_fig,Z_fig,'alphadata',single(Z_fig ~= 0))
    set(gca,'color','r')
    axis xy  equal tight;
    colormap gray;
    hold on;
    imagesc(X_fig,Y_fig,Ni_fig,'alphadata',single(Ni_fig).*.25)
    

    if isfield(tiles,'coastline')
        i=1;
        for i=1:length(tiles.coastline{1})
            if  isempty(tiles.coastline{1}{i}); continue; end
            plot(tiles.coastline{1}{i}(1,:),...
                tiles.coastline{1}{i}(2,:),'b','linewidth',1)
        end
    end
    
    plot([tiles.x0,tiles.x0,tiles.x1,tiles.x1,tiles.x0], [tiles.y0,tiles.y1,tiles.y1,tiles.y0,tiles.y0],'w','linewidth',2)
    
    set(gca,'xlim',[min(X_fig)-500 max(X_fig)+500],'ylim',[min(Y_fig)-500 max(Y_fig)+500]);
    
    
    %set(gcf,'units','normalized');
    %set(gcf,'position',[0.01,0.01,.35,.9])
    qc=load([fileparts(fileName),'/qc.mat']);
    
    fileNames = qc.fileNames;
    
   % alter paths in database if set
    if ~isempty(changePath)
        fileNames = strrep(fileNames,'/mnt/pgc',changePath);
        fileNames = strrep(fileNames,'/','\');        
    end
    
    [~,IA]=intersect(fileNames, fileName);
    
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
        
        if exist('I_dem', 'var')
            clear I_dem;
        end
        
        if flag == 3
            fprintf('entering manual edit mode\n')
            if ~exist(rawfileName, 'file')
                fprintf('** filter mode is NOT available for this image **\n')
            end
            
            j=1;
            while j
                
%                 [~,qc.x{IA}{j},qc.y{IA}{j}] =  roipoly;
                [roi_BW,roi_x,roi_y] =  roipoly;
                polys = {[roi_x roi_y]};
                
                if ~isempty(roi_BW)
                    poly_group = plot(roi_x,roi_y,'g','linewidth',2);

                    if exist(rawfileName, 'file')
                        while j
                            s=input('enter filter mode? (y/n)\n','s');

                            if ~strcmpi(s,'y') && ~strcmpi(s,'n')
                                fprintf('%s not recogized, try again\n',s);
                            else
                                break
                            end
                        end

                        if strcmpi(s,'y')
                            delete(poly_group);
                            poly_group = plot(roi_x,roi_y,'y','linewidth',2);
                            polys = {};
                            
                            % crop figure-size ROI mask to DEM extent
                            roi_BW = roi_BW(Z_r0:Z_r1, Z_c0:Z_c1);
                            
                            % read raw DEM
                            if ~exist('I_dem', 'var')
                                I_dem = readGeotiff(rawfileName);
                            end
                            F = guihandles(fig);
                            F.z = I_dem.z;
                            F.x = I_dem.x;
                            F.y = I_dem.y;
                            
                            % crop ROI and DEM to ROI extent
                            [roi_BW,rcrop,ccrop] = cropzeros(roi_BW);
                            F.z = F.z(rcrop(1):rcrop(2), ccrop(1):ccrop(2));
                            F.x = F.x(ccrop(1):ccrop(2));
                            F.y = F.y(rcrop(1):rcrop(2));
                            F.roi = roi_BW;
                            F.nodata = isnan(F.z);
                            
                            bgcolor = fig.Color;
                            
                            
                            % make slider for kernel size
                            F.kern_val = 5;
                            F.kern_title = uicontrol('Parent',fig,'Style','text','Position',[260,25,100,23],...
                                            'String','Kernel Size','BackgroundColor',bgcolor);
                            F.kern_slide = uicontrol('Parent',fig,'Style','slider','Position',[101,54,419,23],...
                                            'Value',F.kern_val, 'min',1, 'max',20, 'SliderStep',1);
                            F.kern_lmin = uicontrol('Parent',fig,'Style','text','Position',[50,50,43,23],...
                                            'String',1,'BackgroundColor',bgcolor,...
                                            'HorizontalAlignment','Right');
                            F.kern_lmax = uicontrol('Parent',fig,'Style','text','Position',[528,50,43,23],...
                                            'String',[num2str(round(kern_max)),' m'],'BackgroundColor',bgcolor,...
                                            'HorizontalAlignment','Left');
                            F.kern_slide.Callback = @(es,ed) callback_kern_slide(fig, es.Value);
                            
                            
%                             % make slider for slope filter
%                             F.group_slope = hggroup;
%                             [dx,dy] = gradient(z,x,y);
%                             [~,F.slope_raw] = cart2pol(dx,dy);
%                             clear(dx, dy);
                            
                            
                            % make slider for elevation filter
                            F.group_elev = hggroup;
                            F.elev_raw = F.z;
                            F.elev_on = 0;
                            
                            roi_elev = F.elev_raw;
                            roi_elev(~F.roi) = NaN;
                            elev_min = min(roi_elev(:));
                            elev_max = max(roi_elev(:));
                            F.elev_val = (elev_min+elev_max)/2;
                            
                            F.elev_title = uicontrol('Parent',fig,'Style','text','Position',[260,25,100,23],...
                                            'String','Elevation Filter','BackgroundColor',bgcolor);
                            F.elev_togg = uicontrol('Parent',fig,'Style','togglebutton','Position',[20,54,23,23],...
                                            'Value',F.elev_on, 'min',0, 'max',1);
                            F.elev_slide = uicontrol('Parent',fig,'Style','slider','Position',[101,54,419,23],...
                                            'Value',F.elev_val, 'min',elev_min, 'max',elev_max);
                            F.elev_lmin = uicontrol('Parent',fig,'Style','text','Position',[50,50,43,23],...
                                            'String',[num2str(round(elev_min)),' m'],'BackgroundColor',bgcolor,...
                                            'HorizontalAlignment','Right');
                            F.elev_lmax = uicontrol('Parent',fig,'Style','text','Position',[528,50,43,23],...
                                            'String',[num2str(round(elev_max)),' m'],'BackgroundColor',bgcolor,...
                                            'HorizontalAlignment','Left');
                            F.elev_slide.Callback = @(es,ed) callback_elev_slide(fig, es.Value);
                            F.elev_togg.Callback = @(es,ed) callback_elev_togg(fig, es.Value);
                            
                            F.mask = false(size(F.z));
                            prompt = 'apply filters?';
                            while j
                                F.mask_elev = false(size(F.z));
                                F.mask_slope = false(size(F.z));
                                
                                guidata(fig, F) 
                                while j
                                    s=input([prompt,' (y/n)\n'],'s');

                                    if ~strcmpi(s,'y') && ~strcmpi(s,'n')
                                        fprintf('%s not recogized, try again\n',s);
                                    else
                                        break
                                    end
                                end
                                F = guidata(fig);

                                if strcmpi(s,'y')
                                    mask_additions = true(size(F.z));
                                    if F.elev_on
                                        mask_additions = (mask_additions & F.mask_elev);
                                    end
                                    if F.slope_on
                                        mask_additions = (mask_additions & F.mask_slope);
                                    end
                                    F.mask = F.roi & (F.mask | mask_additions);
                                    
                                    prompt = 'apply more filters?';
                                else
                                    B = bwboundaries(F.mask, 'noholes');
                                    polys = [polys; B];
                                    break
                                end
                            end
                            
                            delete(F.group_elev);
                            
                            delete(poly_group);
                            poly_group = hggroup;
                            for k = 1:length(polys)
                               boundary = polys{k};
                               plot(F.x(boundary(:,2)), F.y(boundary(:,1)), 'g','linewidth',2,'Parent',poly_group)
                            end
                        end
                    end
                    
                    if exist('F', 'var')
                        fields = fieldnames(F);
                        for k = 1:numel(fields)
                            delete(F.(fields{i}));
                        end
                        clear F;
                    end
                    
                    if ~isempty(polys)
                        while j
                            s=input('apply mask? (y/n)\n','s');

                            if ~strcmpi(s,'y') && ~strcmpi(s,'n')
                                fprintf('%s not recogized, try again\n',s);
                            else
                                break
                            end
                        end

                        if strcmpi(s,'y')
                            for poly_num = 1:length(polys)
                                qc.x{IA}{j} = polys{poly_num}(:,1);
                                qc.y{IA}{j} = polys{poly_num}(:,2);
                                j = j + 1;
                            end
                        else
                            children = get(poly_group, 'children');
                            if ~isempty(children)
                                delete(children(:));
                            end
                            delete(poly_group);
                        end
                    end
                end
                
                while j
                    s=input('continue editing this image? (y/n)\n','s');
                    
                    if ~strcmpi(s,'y') && ~strcmpi(s,'n')
                        fprintf('%s not recogized, try again\n',s);
                    else
                        break
                    end
                end
                
                if strcmpi(s,'n'); break; end
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
end


function callback_elev_togg(fig, val)
    F = guidata(fig);
    F.elev_on = val;
    guidata(fig, F);
    callback_elev_slide(fig, F.elev_val)
end


function callback_elev_slide(fig, val)
    F = guidata(fig);
    
    children = get(F.group_elev, 'children');
    if ~isempty(children)
        delete(children(:));
    end
    
    if F.elev_on
        if val ~= F.elev_val
            F.elev_val = val;
            F.mask_elev = F.roi & ((F.elev_raw < F.elev_val) | F.nodata);
        end
        B = bwboundaries(F.mask_elev, 'noholes');
        for k = 1:length(B)
           boundary = B{k};
           plot(F.x(boundary(:,2)), F.y(boundary(:,1)), 'cyan','linewidth',2,'Parent',F.group_elev)
        end
    end
    
    guidata(fig, F);
end


function [r] = grade(z,x,y)
    [dx,dy] = gradient(z,x,y);
    [~,r] = cart2pol(dx,dy);
end


function [A,r,c] = cropzeros(varargin)
% cropnans crop array of bordering nans
%
% [A,r,c] = cropnans(A)

A=varargin{1};
buff=0;
if nargin == 2; buff=varargin{2}; end

r = [];
c = [];


% M = ~isnan(A);

if ~any(A(:)); return; end

rowsum = sum(A) ~= 0;
colsum = sum(A,2) ~= 0;

c(1) = find(rowsum,1,'first')-buff;
c(2) = find(rowsum,1,'last')+buff;

r(1) = find(colsum,1,'first')-buff;
r(2) = find(colsum,1,'last')+buff;

if c(1) < 1; c(1)=1; end
if r(1) < 1; r(1)=1; end
if c(2) > size(A,2); c(2)=size(A,2); end
if r(2) > size(A,1); r(2)=size(A,1); end

A = A(r(1):r(2),c(1):c(2));
end
