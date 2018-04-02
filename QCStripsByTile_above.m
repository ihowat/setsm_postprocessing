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
    
    
    image_fig = figure(1);
    image_ax = gca;
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
                [roi_BW,roi_x,roi_y] = roipoly;
                
                if ~isempty(roi_BW)
                    polys = {[roi_x roi_y]};
                    poly_group = plot(roi_x,roi_y,'g','linewidth',2);

                    if exist(rawfileName, 'file')
                        while j
                            s=input('enter filter mode? (1/0)\n','s');
                            if ~strcmpi(s,'1') && ~strcmpi(s,'0')
                                fprintf('%s not recogized, try again\n',s);
                            else
                                break
                            end
                        end

                        if strcmpi(s,'1')
                            % Read raw DEM.
                            if ~exist('I_dem', 'var')
                                I_dem = readGeotiff(rawfileName);
                            end
                            % Crop figure-size ROI mask to DEM extent.
                            roi_BW = roi_BW(Z_r0:Z_r1, Z_c0:Z_c1);
                            % Clear ROI poly in prep for filter mode.
                            delete(poly_group);
                            % Enter filter mode.
                            [polys, poly_group] = filter_mode(image_fig, image_ax, I_dem, roi_BW, roi_x, roi_y);
                        end
                    end
                    
                    if ~isempty(polys)
                        while j
                            s=input('apply mask? (1/0)\n','s');
                            if ~strcmpi(s,'1') && ~strcmpi(s,'0')
                                fprintf('%s not recogized, try again\n',s);
                            else
                                break
                            end
                        end

                        if strcmpi(s,'1')
                            % Write poly(s) to qc file.
                            for poly_num = 1:length(polys)
                                qc.x{IA}{j} = polys{poly_num}(:,1);
                                qc.y{IA}{j} = polys{poly_num}(:,2);
                                j = j + 1;
                            end
                        else
                            % Clear polys created with this edit.
                            children = get(poly_group, 'children');
                            if ~isempty(children)
                                delete(children(:));
                            end
                            delete(poly_group);
                        end
                    end
                end
                
                while j
                    s=input('continue editing this image? (1/0)\n','s');
                    
                    if ~strcmpi(s,'1') && ~strcmpi(s,'0')
                        fprintf('%s not recogized, try again\n',s);
                    else
                        break
                    end
                end
                
                if strcmpi(s,'0'); break; end
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


function [A, entropy_kernel_size] = get_pixvals_raw(F, item_num)
    entropy_kernel_size = [];
    switch item_num
        case F.nans.item_num
            A = isnan(F.z);
        case F.slope.item_num
            [dx,dy] = gradient(F.z,F.x,F.y);
            [~,A] = cart2pol(dx,dy);
            A = A*100;  % get in terms of percent grade
        case F.entropy.item_num
            entropy_kernel_size = F.kernel_smooth.ui_slider.Value;
            if mod(entropy_kernel_size, 2) == 0
                entropy_kernel_size = entropy_kernel_size + 1;
            end
            z_uint16 = cast(F.z, 'uint16');
            z_subtraction =  movmax(z_uint16, entropy_kernel_size) - movmin(z_uint16, entropy_kernel_size);
            A = entropyfilt(z_subtraction, true(entropy_kernel_size));
        case F.elev.item_num
            A = F.z;
        otherwise
            error('`get_pixvals_raw` item_num=%d not recognized\n', item_num);
    end
end


function [E, item_name] = get_struct_element_by_item_num(F, item_num)
    struct_fields = fieldnames(F);
    for i = 1:numel(struct_fields)
        E = F.(struct_fields{i});
        if isstruct(E)
            element_fields = fieldnames(E);
            for j = 1:numel(element_fields)
                if strcmp(element_fields{j}, 'item_num')
                    if E.item_num == item_num
                        item_name = struct_fields{i};
                        return
                    end
                end
            end
        end
    end
    error('`get_struct_element_by_item_num` unsuccessfull` for item_num=%d', item_num);
end


function [polys, poly_group] = filter_mode(fig, ax, I_dem, roi_BW, roi_x, roi_y)
    polys = {};
    poly_group = plot(roi_x,roi_y,'y','linewidth',2);

    F = guihandles(fig);
    F.z = I_dem.z;
    F.x = I_dem.x;
    F.y = I_dem.y;

    % Crop ROI and DEM to ROI extent.
    [roi_BW,rcrop,ccrop] = cropzeros(roi_BW);
    F.z = F.z(rcrop(1):rcrop(2), ccrop(1):ccrop(2));
    F.x = F.x(ccrop(1):ccrop(2));
    F.y = F.y(rcrop(1):rcrop(2));
    F.roi = roi_BW;
    F.basemask_array = false(size(F.z));

    
    %%% UI CONTROL PANEL %%%
    fig_panel = uifigure;
    
    %% NANS %%
    E = struct();
    E.item_num = 1;
    E.item_group = 1;
    E.plot_group = hggroup;
    E.plot_color = 'magenta';
    E.mask_kernel_size = 1;
    E.ui_title = uilabel(fig_panel, 'Position',[155 362 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','Mask NaNs');
    E.ui_toggle = uibutton(fig_panel, 'state', 'Position',[115 355 35 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_on(fig, E.item_num, btn, false),...
        'Value',false, 'Text','OFF');  % DEFAULT
    E.ui_display = uibutton(fig_panel, 'state', 'Position',[430 355 50 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_display(fig, E.item_num, btn),...
        'Value',true, 'Text','SHOW');  % DEFAULT
    F.nans = E;
    
    %% SMOOTH KERNEL SIZE %%
    E = struct();
    E.item_num = 2;
    E.item_group = 1;
    max_kernel_size = 25;  % SETTING
    E.ui_title = uilabel(fig_panel, 'Position',[155 300 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','Smooth Kernel Size (sq. side length)');
    E.ui_toggle = uibutton(fig_panel, 'state', 'Position',[115 275 35 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_on(fig, E.item_num, btn, true),...
        'Value',false, 'Text','OFF');  % DEFAULT
    E.ui_slider_dec = uibutton(fig_panel, 'push', 'Position',[160 280 20 25], 'Text','-',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, -1, true, true));
    E.ui_slider = uislider(fig_panel, 'Position',[190 290 200 3],...
        'ValueChangedFcn',@(sld,event) slider_drag_integer(fig, E.item_num),...
        'MajorTicksMode','manual', 'MajorTicks',[0:5:max_kernel_size],...
        'MinorTicksMode','manual', 'MinorTicks',[1:1:max_kernel_size],...
        'Limits',[1 max_kernel_size], 'Value',5);  % DEFAULT
    E.ui_slider_inc = uibutton(fig_panel, 'push', 'Position',[400 280 20 25], 'Text','+',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, 1, true, true));
    F.kernel_smooth = E;
    
    %% SLOPE FILTER %%
    E = struct();
    E.item_num = 3;
    E.item_group = 1;
    E.plot_group = hggroup;
    E.plot_color = 'green';
    E.mask_kernel_size = 1;
    E.ui_title = uilabel(fig_panel, 'Position',[155 220 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','Slope Filter (% grade)');
    E.ui_logic = uibuttongroup(fig_panel, 'Position',[20 192 50 40],...
        'SelectionChangedFcn',@(bg,event) recompute_basemask(fig, E.item_num));
    E.ui_logic_or = uiradiobutton(E.ui_logic, 'Position',[3 7 45 30], 'Text','OR');
    E.ui_logic_and = uiradiobutton(E.ui_logic, 'Position',[3 -11 45 30], 'Text','AND');
    E.ui_logic.SelectedObject = E.ui_logic_or;  % DEFAULT
    E.mask_logic = E.ui_logic.SelectedObject;
    E.ui_cond = uibuttongroup(fig_panel, 'Position',[70 192 35 40],...
        'SelectionChangedFcn',@(bg,event) recompute_basemask(fig, E.item_num));
    E.ui_cond_gt = uiradiobutton(E.ui_cond, 'Position',[5 7 30 30], 'Text','>');
    E.ui_cond_lt = uiradiobutton(E.ui_cond, 'Position',[5 -11 30 30], 'Text','<');
    E.ui_cond.SelectedObject = E.ui_cond_gt;  % DEFAULT
    E.mask_cond = E.ui_cond.SelectedObject;
    E.ui_toggle = uibutton(fig_panel, 'state', 'Position',[115 195 35 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_on(fig, E.item_num, btn, false),...
        'Value',false, 'Text','OFF');  % DEFAULT
    E.ui_slider_dec = uibutton(fig_panel, 'push', 'Position',[160 200 20 25], 'Text','-',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, -1, false));
    E.ui_slider = uislider(fig_panel, 'Position',[190 210 200 3],...
        'ValueChangedFcn',@(sld,event) recompute_basemask(fig, E.item_num),...
        'Limits',[0 100], 'Value',0);  % DEFAULT
    E.ui_slider_inc = uibutton(fig_panel, 'push', 'Position',[400 200 20 25], 'Text','+',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, 1, false));
    E.mask_slider_value = E.ui_slider.Value;
    E.ui_display = uibutton(fig_panel, 'state', 'Position',[430 195 50 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_display(fig, E.item_num, btn),...
        'Value',true, 'Text','SHOW');  % DEFAULT
    F.slope = E;
    
    %% ENTROPY FILTER %%
    E = struct();
    E.item_num = 4;
    E.item_group = 1;
    E.plot_group = hggroup;
    E.plot_color = 'green';
    E.mask_kernel_size = 1;
    E.ui_title = uilabel(fig_panel, 'Position',[155 140 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','Entropy Filter');
    E.ui_logic = uibuttongroup(fig_panel, 'Position',[20 112 50 40],...
        'SelectionChangedFcn',@(bg,event) recompute_basemask(fig, E.item_num));
    E.ui_logic_or = uiradiobutton(E.ui_logic, 'Position',[3 7 45 30], 'Text','OR');
    E.ui_logic_and = uiradiobutton(E.ui_logic, 'Position',[3 -11 45 30], 'Text','AND');
    E.ui_logic.SelectedObject = E.ui_logic_or;  % DEFAULT
    E.mask_logic = E.ui_logic.SelectedObject;
    E.ui_cond = uibuttongroup(fig_panel, 'Position',[70 112 35 40],...
        'SelectionChangedFcn',@(bg,event) recompute_basemask(fig, E.item_num));
    E.ui_cond_gt = uiradiobutton(E.ui_cond, 'Position',[5 7 30 30], 'Text','>');
    E.ui_cond_lt = uiradiobutton(E.ui_cond, 'Position',[5 -11 30 30], 'Text','<');
    E.ui_cond.SelectedObject = E.ui_cond_gt;  % DEFAULT
    E.mask_cond = E.ui_cond.SelectedObject;
    E.ui_toggle = uibutton(fig_panel, 'state', 'Position',[115 115 35 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_on(fig, E.item_num, btn, false),...
        'Value',false, 'Text','OFF');  % DEFAULT
    E.ui_slider_dec = uibutton(fig_panel, 'push', 'Position',[160 120 20 25], 'Text','-',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, -1, false));
    E.ui_slider = uislider(fig_panel, 'Position',[190 130 200 3],...
        'ValueChangedFcn',@(sld,event) recompute_basemask(fig, E.item_num),...
        'Limits',[0 10], 'Value',0);  % DEFAULT
    E.ui_slider_inc = uibutton(fig_panel, 'push', 'Position',[400 120 20 25], 'Text','+',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, 1, false));
    E.mask_slider_value = E.ui_slider.Value;
    E.ui_display = uibutton(fig_panel, 'state', 'Position',[430 115 50 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_display(fig, E.item_num, btn),...
        'Value',true, 'Text','SHOW');  % DEFAULT
    F.entropy = E;
    
    %% ELEVATION FILTER %%
    E = struct();
    E.item_num = 5;
    E.item_group = 1;
    E.plot_group = hggroup;
    E.plot_color = 'cyan';
    E.mask_kernel_size = 1;
    E.ui_title = uilabel(fig_panel, 'Position',[155 60 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','Elevation Filter (m)');
    E.ui_logic = uibuttongroup(fig_panel, 'Position',[20 32 50 40],...
        'SelectionChangedFcn',@(bg,event) recompute_basemask(fig, E.item_num));
    E.ui_logic_or = uiradiobutton(E.ui_logic, 'Position',[3 7 45 30], 'Text','OR');
    E.ui_logic_and = uiradiobutton(E.ui_logic, 'Position',[3 -11 45 30], 'Text','AND');
    E.ui_logic.SelectedObject = E.ui_logic_or;  % DEFAULT
    E.mask_logic = E.ui_logic.SelectedObject;
    E.ui_cond = uibuttongroup(fig_panel, 'Position',[70 32 35 40],...
        'SelectionChangedFcn',@(bg,event) recompute_basemask(fig, E.item_num));
    E.ui_cond_gt = uiradiobutton(E.ui_cond, 'Position',[5 7 30 30], 'Text','>');
    E.ui_cond_lt = uiradiobutton(E.ui_cond, 'Position',[5 -11 30 30], 'Text','<');
    E.ui_cond.SelectedObject = E.ui_cond_gt;  % DEFAULT
    E.mask_cond = E.ui_cond.SelectedObject;
    E.ui_toggle = uibutton(fig_panel, 'state', 'Position',[115 35 35 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_on(fig, E.item_num, btn, false),...
        'Value',false, 'Text','OFF');  % DEFAULT
    E.ui_slider_dec = uibutton(fig_panel, 'push', 'Position',[160 40 20 25], 'Text','-',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, -1, false));
    E.ui_slider = uislider(fig_panel, 'Position',[190 50 200 3],...
        'ValueChangedFcn',@(sld,event) recompute_basemask(fig, E.item_num),...
        'Limits',[0 100], 'Value',0);  % DEFAULT
    E.ui_slider_inc = uibutton(fig_panel, 'push', 'Position',[400 40 20 25], 'Text','+',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, 1, false));
    E.mask_slider_value = E.ui_slider.Value;
    E.ui_display = uibutton(fig_panel, 'state', 'Position',[430 35 50 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_display(fig, E.item_num, btn),...
        'Value',true, 'Text','SHOW');  % DEFAULT
    F.elev = E;
    
    
    %% RESULT MASK %%
    E = struct();
    E.item_num = 6;
    E.item_group = 2;
    E.mask_array = false(size(F.basemask_array));
    E.plot_group = hggroup;
    E.plot_color = 'cyan';
    E.ui_title = uilabel(fig_panel, 'Position',[655 362 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','RESULT MASK');
    E.ui_cond = uibuttongroup(fig_panel, 'Position',[570 352 35 40],...
        'SelectionChangedFcn',@(bg,event) recompute_basemask(fig, E.item_num));
    E.ui_cond_gt = uiradiobutton(E.ui_cond, 'Position',[5 7 30 30], 'Text','>');
    E.ui_cond_lt = uiradiobutton(E.ui_cond, 'Position',[5 -11 30 30], 'Text','<');
    E.ui_cond.SelectedObject = E.ui_cond_gt;  % DEFAULT
    E.mask_cond = E.ui_cond.SelectedObject;
    E.ui_display = uibutton(fig_panel, 'state', 'Position',[930 355 50 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_display(fig, E.item_num, btn),...
        'Value',true, 'Text','SHOW');  % DEFAULT
    F.result = E;
    
    %% CLUSTER SIZE FILTER %%
    E = struct();
    E.item_num = 7;
    E.item_group = 2;
    E.plot_group = hggroup;
    E.plot_color = 'cyan';
    E.ui_title = uilabel(fig_panel, 'Position',[655 300 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','Cluster Size Filter (px)');
    E.ui_cond = uibuttongroup(fig_panel, 'Position',[570 272 35 40],...
        'SelectionChangedFcn',@(bg,event) recompute_basemask(fig, E.item_num));
    E.ui_cond_gt = uiradiobutton(E.ui_cond, 'Position',[5 7 30 30], 'Text','>');
    E.ui_cond_lt = uiradiobutton(E.ui_cond, 'Position',[5 -11 30 30], 'Text','<');
    E.ui_cond.SelectedObject = E.ui_cond_gt;  % DEFAULT
    E.mask_cond = E.ui_cond.SelectedObject;
    E.ui_toggle = uibutton(fig_panel, 'state', 'Position',[615 275 35 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_on(fig, E.item_num, btn, false),...
        'Value',false, 'Text','OFF');  % DEFAULT
    E.ui_slider_dec = uibutton(fig_panel, 'push', 'Position',[660 280 20 25], 'Text','-',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, -1, true));
    E.ui_slider = uislider(fig_panel, 'Position',[690 290 200 3],...
        'ValueChangedFcn',@(sld,event) recompute_resultmask(fig, E.item_num, flase, false, true),...
        'Limits',[0 100], 'Value',0);  % DEFAULT
    E.ui_slider_inc = uibutton(fig_panel, 'push', 'Position',[900 280 20 25], 'Text','+',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, 1, true));
    E.mask_slider_value = E.ui_slider.Value;
    F.cluster_size = E;
    
    %% CLUSTER MERGE DISTANCE %%
    E = struct();
    E.item_num = 8;
    E.item_group = 2;
    E.plot_group = hggroup;
    E.plot_color = 'cyan';
    E.ui_title = uilabel(fig_panel, 'Position',[655 220 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','Cluster Merge Distance (m)');
    E.ui_toggle = uibutton(fig_panel, 'state', 'Position',[615 195 35 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_on(fig, E.item_num, btn, false),...
        'Value',false, 'Text','OFF');  % DEFAULT
    E.ui_slider_dec = uibutton(fig_panel, 'push', 'Position',[660 200 20 25], 'Text','-',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, -1, false));
    E.ui_slider = uislider(fig_panel, 'Position',[690 210 200 3],...
        'ValueChangedFcn',@(sld,event) recompute_basemask(fig, E.item_num),...
        'Limits',[0 100], 'Value',0);  % DEFAULT
    E.ui_slider_inc = uibutton(fig_panel, 'push', 'Position',[900 200 20 25], 'Text','+',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, 1, false));
    E.mask_slider_value = E.ui_slider.Value;
    F.cluster_dist = E;
    
    %% CLUSTER CONCAVITY %%
    E = struct();
    E.item_num = 9;
    E.item_group = 2;
    E.plot_group = hggroup;
    E.plot_color = 'cyan';
    E.ui_title = uilabel(fig_panel, 'Position',[655 140 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','Cluster Concavity');
    E.ui_toggle = uibutton(fig_panel, 'state', 'Position',[615 115 35 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_on(fig, E.item_num, btn, false),...
        'Value',false, 'Text','OFF');  % DEFAULT
    E.ui_slider_dec = uibutton(fig_panel, 'push', 'Position',[660 120 20 25], 'Text','-',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, -1, false));
    E.ui_slider = uislider(fig_panel, 'Position',[690 130 200 3],...
        'ValueChangedFcn',@(sld,event) recompute_basemask(fig, E.item_num),...
        'Limits',[0 1], 'Value',0);  % DEFAULT
    E.ui_slider_incv = uibutton(fig_panel, 'push', 'Position',[900 120 20 25], 'Text','+',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, 1, false));
    E.mask_slider_value = E.ui_slider.Value;
    F.cluster_concav = E;
    
    %% DILATION SIZE %%
    E = struct();
    E.item_num = 10;
    E.item_group = 2;
    max_kernel_size = 25;  % SETTING
    E.ui_title = uilabel(fig_panel, 'Position',[655 60 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','Dilation Kernel Size (sq. side length)');
    E.ui_toggle = uibutton(fig_panel, 'state', 'Position',[615 35 35 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_on(fig, E.item_num, btn, true),...
        'Value',false, 'Text','OFF');  % DEFAULT
    E.ui_slider_dec = uibutton(fig_panel, 'push', 'Position',[660 40 20 25], 'Text','-',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, -1, true, true));
    E.ui_slider = uislider(fig_panel, 'Position',[690 50 200 3],...
        'ValueChangedFcn',@(sld,event) slider_drag_integer(fig, E.item_num),...
        'MajorTicksMode','manual', 'MajorTicks',[0:5:max_kernel_size],...
        'MinorTicksMode','manual', 'MinorTicks',[1:1:max_kernel_size],...
        'Limits',[1 max_kernel_size], 'Value',5);  % DEFAULT
    E.ui_slider_inc = uibutton(fig_panel, 'push', 'Position',[900 40 20 25], 'Text','+',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, 1, true, true));
    F.kernel_dilate = E;
    
    
    guidata(fig, F);
    
    prompt = 'apply filters?';
    while true
        guidata(fig, F);
        while true
            s=input([prompt,' (1/0)\n'],'s');
            if ~strcmpi(s,'1') && ~strcmpi(s,'0')
                fprintf('%s not recogized, try again\n',s);
            else
                break
            end
        end
        F = guidata(fig);
        
%         if strcmpi(s,'1')
%         end
    end
end


function toggle_item_on(fig, item_num, btn, isKernel)
    if ~exist('isKernel', 'var') || isempty(isKernel)
            isKernel = false;
    end
    
    F = guidata(fig);
    [E,item_name] = get_struct_element_by_item_num(F, item_num);
    
    if btn.Value == 0
        btn.Text = 'OFF';
    else
        btn.Text = 'ON';
    end
    guidata(fig, F);
    
    if E.item_group == 1
        if isKernel
            recompute_basemask(fig, item_num);
        else
            recompute_basemask(fig, item_num, true, true);
        end
    else
        recompute_resultmask(fig, item_num, false, true);
    end
end


function toggle_item_display(fig, item_num, btn)
    F = guidata(fig);
    if btn.Value == 0
        btn.Text = 'HIDE';
    else
        btn.Text = 'SHOW';
    end
    guidata(fig, F);
    redraw_item(fig, item_num, false);
end


function slider_step(fig, item_num, direction, roundToInt, isKernel)
    if ~exist('roundToInt', 'var') || isempty(roundToInt)
            roundToInt = false;
    end
    if ~exist('isKernel', 'var') || isempty(isKernel)
            isKernel = false;
    end
    
    F = guidata(fig);
    [E,item_name] = get_struct_element_by_item_num(F, item_num);
    
    slider_limits = E.ui_slider.Limits;
    slider_min = slider_limits(1,1);
    slider_max = slider_limits(1,2);
    
    if isKernel
        step_size = 1;
    else
        step_size = (slider_max-slider_min)/20;
    end
    
    old_slider_value = E.ui_slider.Value;
    new_slider_value = min(slider_max, max(slider_min, old_slider_value + direction*step_size));
    if roundToInt
        new_slider_value = round(new_slider_value);
    end
    if new_slider_value == old_slider_value
        return;
    end
    
    E.ui_slider.Value = new_slider_value;
    
    F.(item_name) = E;
    guidata(fig, F);
    
    if E.item_group == 1
        recompute_basemask(fig, item_num);
    else
        recompute_resultmask(fig, item_num, false, false, true);
    end
end

    
function slider_drag_integer(fig, item_num)
    F = guidata(fig);
    [E,item_name] = get_struct_element_by_item_num(F, item_num);
    
    E.ui_slider.Value = round(E.ui_slider.Value);
    
    F.(item_name) = E;
    guidata(fig, F);
    
    if E.item_group == 1
        recompute_basemask(fig, item_num);
    else
        recompute_resultmask(fig, item_num, false, false, true);
    end
end


function redraw_item(fig, item_num, modified)
    F = guidata(fig);
    [E,item_name] = get_struct_element_by_item_num(F, item_num);
    
    if ~isfield(E, 'plot_group')
        return;
    end
    
    % Clear drawn polygons for this item.
    prev_plot_polys = get(E.plot_group, 'children');
    if ~isempty(prev_plot_polys)
        delete(prev_plot_polys(:));
    end
    
    if E.ui_toggle.Value == 1 && E.ui_display.Value == 1
        % Draw polygons for this item.
        if modified || ~isfield(E, 'polys')
            E.polys = bwboundaries(E.mask_array, 'noholes');
        end
        for k = 1:length(E.polys)
           poly = E.polys{k};
           plot(F.x(poly(:,2)), F.y(poly(:,1)), E.plot_color,'linewidth',2, 'Parent',E.plot_group)
        end
    end
    
    F.(item_name) = E;
    guidata(fig, F);
end


function [modified] = cluster_filter_size(basemask_modified, slider_modified)
    F = guidata(fig);
    E = F.cluster_size;
    
    modified = basemask_modified | slider_modified;
    
    latest_prefilt_array = F.basemask_array;
    
    if ~isfield(E, 'prefilt_array') || ~isfield(F.basemask_RP, 'Area')
        basemask_modified = true;
        modified = true;
    end
    if E.mask_cond ~= E.ui_cond.SelectedObject
        E.mask_cond = E.ui_cond.SelectedObject;
        modified = true;
    end
    if E.mask_slider_value ~= E.ui_slider.Value
        E.mask_slider_value = E.ui_slider.Value;
        modified = true;
    end
    if ~basemask_modified && ~isequal(latest_prefilt_array, E.prefilt_array)
        basemask_modified = true;
        modified = true;
    end
    if modified
        if basemask_modified
            E.prefilt_array = F.basemask_array;
            rp = regionprops(F.basemask_CC, 'Area');
            F.basemask_RP.Area = [rp.Area];
            [new_slider_min,~] = min(F.basemask_RP.Area);
            [new_slider_max,~] = max(F.basemask_RP.Area);
            if new_slider_min == new_slider_max
                new_slider_max = new_slider_min + 1;
            end
            new_slider_limits = zeros(size(1,2));
            new_slider_limits(1,1) = new_slider_min;
            new_slider_limits(1,2) = new_slider_max;
            E.ui_slider.Limits = new_slider_limits;
        end
        if E.ui_cond.SelectedObject == E.ui_cond_lt
            new_postfilt_CC_ind = find(F.basemask_RP.Area < E.ui_slider.Value);
        else
            new_postfilt_CC_ind = find(F.basemask_RP.Area > E.ui_slider.Value);
        end
        if isfield(E, 'postfilt_CC_ind') && slider_modified && ~basemask_modified && isequal(new_postfilt_CC_ind, E.postfilt_CC_ind)
            modified = false;
        else
            E.postfilt_CC_ind = new_postfilt_CC_ind;
        end
        if modified
            new_postfilt_array = bwareaopen(E.prefilt_array, E.ui_slider.Value);
            if E.ui_cond.SelectedObject == E.ui_cond_lt
                new_postfilt_array = xor(E.prefilt_array, new_postfilt_array);
            end
            E.postfilt_array = new_postfilt_array;
        end
        F.cluster_size = E;
        guidata(fig, F);
    end
end


function cluster_filter_dist(basemask_modified, slider_modified, dilated)
    F = guidata(fig);
    E = F.cluster_dist;
    
    modified = basemask_modified | slider_modified | dilated;
    sizefilt_modified = false;
    
    if F.kernel_dilate.ui_toggle.Value == 1
        latest_prefilt_array = F.kernel_dilate.postfilt_array;
    elseif F.cluster_size.ui_toggle.Value == 1
        latest_prefilt_array = F.cluster_size.postfilt_array;
    else
        latest_prefilt_array = F.basemask_array;
    end
    
    if F.cluster_size.ui_toggle.Value == 1
        latest_prefilt_CC_ind = F.cluster_size.postfilt_CC_ind;
    else
        latest_prefilt_CC_ind = [1:F.basemask_CC.NumObjects];
    end
    if ~isfield(E, 'prefilt_CC_ind') || ~isequal(E.prefilt_CC_ind, latest_prefilt_CC_ind)
        E.prefilt_CC_ind = latest_prefilt_CC_ind;
        sizefilt_modified = true;
        modified = true;
    end
    
    if ~isfield(F.basemask_RP, 'cent_coords')
        basemask_modified = true;
        modified = true;
    end
    if E.mask_slider_value ~= E.ui_slider.Value
        E.mask_slider_value = E.ui_slider.Value;
        modified = true;
    end
    if ~basemask_modified && isfield(E, 'prefilt_CC_ind') && ~isequal(E.prefilt_CC_ind, E.prefilt_array)
        basemask_modified = true;
        modified = true;
    end
    if modified
        if basemask_modified
            rp = regionprops(F.basemask_CC, 'Centroid');
            num_clusters = F.basemask_CC.NumObjects;
            [num_rows, num_cols] = size(latest_prefilt_array);
            cent_coords = round(reshape([rp.Centroid], [2,num_clusters]).');
            cent_coords(cent_coords == 0) = 1;
            cent_coords((cent_coords(:,1) > array_size(2)), 1) = num_cols;
            cent_coords((cent_coords(:,2) > array_size(1)), 1) = num_rows;
            F.basemask_RP.cent_coords = cent_coords;
            cent_coords_col = cent_coords(:,1);
            cent_coords_row = cent_coords(:,2);
            F.basemask_RP.dist_array = sqrt( ...
                power(repmat(cent_coords_row, [1,num_clusters]) - repmat(cent_coords_row.', [num_clusters,1]), 2) ...
              + power(repmat(cent_coords_col, [1,num_clusters]) - repmat(cent_coords_col.', [num_clusters,1]), 2));
            sizefilt_modified = true;
        end
        
        if sizefilt_modified
            E.dist_array = F.basemask_RP.dist_array(latest_prefilt_CC_ind, :);
            new_slider_min = min(F.basemask_RP.dist_array);
            new_slider_max = max(F.basemask_RP.dist_array);
            if new_slider_min == new_slider_max
                new_slider_max = new_slider_min + 1;
            end
            new_slider_limits = zeros(size(1,2));
            new_slider_limits(1,1) = new_slider_min;
            new_slider_limits(1,2) = new_slider_max;
            E.ui_slider.Limits = new_slider_limits;
        end
        
        if ~dilated
            num_clusters = F.basemask_CC.NumObjects;
            dist_array = E.dist_array;
            dist_filter_value = E.mask_slider_value;
            prefilt_CC_ind = E.prefilt_CC_ind;
            new_postfilt_CC_ind = [];
            for i=1:num_clusters
                for j=1:(i-1)
                    if dist_array(i,j) > dist_filter_value
                        new_postfilt_CC_ind = [new_postfilt_CC_ind; [prefilt_CC_ind(i) prefilt_CC_ind(j)]];
                    end
                end
            end
            if isfield(E, 'postfilt_CC_ind') && isequal(new_postfilt_CC_ind, E.postfilt_CC_ind)
                if slider_modified && ~basemask_modified
                    modified = false;
                elseif isfield(E, 'prefilt_array') && isequal(latest_prefilt_array, E.prefilt_array)
                    modified = false;
                end
            else
                E.postfilt_CC_ind = new_postfilt_CC_ind;
            end
        end
        
        if modified
            E.prefilt_array = latest_prefilt_array;
            new_postfilt_array = E.prefilt_array;
            array_size = size(new_postfilt_array);
            postfilt_CC_ind = E.postfilt_CC_ind;
            cent_coords = F.basemask_RP.cent_coords;
            for k=1:length(postfilt_CC_ind)
                cluster_pair_ind = postfilt_CC_ind(k,:);
                i = cluster_pair_ind(1);
                j = cluster_pair_ind(2);
                x = [cent_coords(i,1) cent_coords(j,1)];
                y = [cent_coords(i,2) cent_coords(j,2)];
                nPoints = max(abs(diff(x)), abs(diff(y))) + 1;
                rIndices = round(linspace(y(1), y(2), nPoints));
                cIndices = round(linspace(x(1), x(2), nPoints));
                indices = sub2ind(array_size, rIndex, cIndex);
                new_postfilt_array(indices) = true;
            end
            %%%%%%%%%%%%%%%%%%%%
            if isfield(E, 'postfilt_array') && isequal(new_postfilt_array, E.postfilt_array)
                fprintf('F.cluster_dist.postfilt_array should be different here, but is not!');
            end
            %%%%%%%%%%%%%%%%%%%%
            E.postfilt_array = new_postfilt_array;
        end
        
        F.cluster_dist = E;
        guidata(fig, F);
    end
end


function cluster_dilate(basemask_modified, slider_modified)
    F = guidata(fig);
    E = F.kernel_dilate;
    
    modified = basemask_modified | slider_modified;
    
    if F.cluster_size.ui_toggle.Value == 1
        latest_postfilt_array = F.cluster_size.postfilt_array;
    else
        latest_postfilt_array = F.basemask_array;
    end
    
    if ~isfield(E, 'prefilt_array')
        modified = true;
    end
    if E.mask_slider_value ~= E.ui_slider.Value
        E.mask_slider_value = E.ui_slider.Value;
        modified = true;
    end
    if ~modified && ~isequal(E.prefilt_array, latest_postfilt_array)
        modified = true;
    end
    if modified
        E.prefilt_array = latest_postfilt_array;
        E.postfilt_array = imdilate(E.prefilt_array, true(E.ui_slider.Value));
        F.kernel_dilate = E;
        guidata(fig, F);
    end
end


function recompute_resultmask(fig, item_num, modified, redraw, cluster_slider)
    if ~exist('modified', 'var') || isempty(modified)
        modified = false;
    end
    if ~exist('redraw', 'var') || isempty(redraw)
        redraw = false;
    end
    if ~exist('cluster_slider', 'var') || isempty(cluster_slider)
        cluster_slider = false;
    end
    
    F = guidata(fig);
    
    if (   (   F.cluster_size.ui_toggle.Value == 1 ...
            || F.cluster_merge.ui_toggle.Value == 1 ...
            || F.cluster_concav.ui_toggle.Value == 1) ...
        && (~isfield(F, 'basemask_CC') || item_num < F.cluster_size.item_num))
        F.basemask_CC = bwconncomp(F.basemask_array);
        F.basemask_RP = struct();
        guidata(fig, F);
    end
    
    if any(F.basemask_array(:))
        
        if F.cluster_size.ui_toggle.Value == 1 && (item_num == F.cluster_size.item_num || (item_num < F.cluster_size.item_num && modified))
            modified = cluster_filter_size(modified, (cluster_slider && item_num == F.cluster_size.item_num));
            if item_num == F.cluster_size.item_num && redraw
                modified = true;
            end
        end
        if F.kernel_dilate.ui_toggle.Value == 1 && (item_num == F.kernel_dilate.item_num || (item_num < F.kernel_dilate.item_num && modified))
            modified = cluster_dilate(modified, (cluster_slider && item_num == F.kernel_dilate.item_num));
            if item_num == F.kernel_dilate.item_num && redraw
                modified = true;
            end
        end
        if F.cluster_dist.ui_toggle.Value == 1 && (item_num == F.cluster_dist.item_num || (item_num < F.cluster_dist.item_num && modified))
            modified = cluster_filter_dist((modified && item_num < F.cluster_dist.item_num), (cluster_slider && item_num == F.cluster_dist.item_num), item_num == F.kernel_dilate.item_num);
            if item_num == F.cluster_dist.item_num && redraw
                modified = true;
            end
        end
        
        
    end
    
    
    
end


function [modified] = recompute_basemask_component(fig, item_num, redraw)
    if ~exist('redraw', 'var') || isempty(redraw)
        redraw = false;
    end

    F = guidata(fig);
    [E, item_name] = get_struct_element_by_item_num(F, item_num);
    
    modified = false;
    return_modified = false;
    minmax_modified = false;
    
    if E.ui_toggle.Value == 1
        smooth_entropy = false;
        if (   ~isfield(E, 'pixvals_raw') ...
            || (E.item_num == F.entropy.item_num && redraw && E.entropy_kernel_size ~= F.kernel_smooth.ui_slider.Value))
            [E.pixvals_raw, entropy_kernel_size] = get_pixvals_raw(F, item_num);
            E.pixvals_smooth = E.pixvals_raw;
            minmax_modified = true;
            modified = true;
            if E.item_num == F.entropy.item_num
                E.entropy_kernel_size = entropy_kernel_size;
                E.ui_title.Text = sprintf('Entropy Filter (%dx%d kernel)', entropy_kernel_size, entropy_kernel_size);
                smooth_entropy = true;
            end
        end
        if (F.kernel_smooth.ui_toggle.Value == 1 && E.mask_kernel_size ~= F.kernel_smooth.ui_slider.Value) || smooth_entropy
            if F.kernel_smooth.ui_slider.Value > 1
                if E.item_num ~= F.nans.item_num
                    E.pixvals_smooth = conv2(E.pixvals_raw, ones(F.kernel_smooth.ui_slider.Value)/(F.kernel_smooth.ui_slider.Value.^2), 'same');
                else
                    E.pixvals_smooth = imdilate(E.pixvals_raw, true(F.kernel_smooth.ui_slider.Value));
                end
            else
                E.pixvals_smooth = E.pixvals_raw;
            end
            E.mask_kernel_size = F.kernel_smooth.ui_slider.Value;
            minmax_modified = true;
            modified = true;
        elseif F.kernel_smooth.ui_toggle.Value == 0 && E.mask_kernel_size ~= 1
            E.pixvals_smooth = E.pixvals_raw;
            E.mask_kernel_size = 1;
            minmax_modified = true;
            modified = true;
        end
        if E.item_num ~= F.nans.item_num
            if minmax_modified
                new_slider_min = round(min(E.pixvals_smooth(:)), 2);
                new_slider_max = round(max(E.pixvals_smooth(:)), 2);
                if new_slider_min == new_slider_max
                    new_slider_max = new_slider_min + 0.01;
                end
%                 new_range = new_slider_max - new_slider_min;
%                 new_slider_min = new_slider_min - round(new_range*0.05, 3);
%                 new_slider_max = new_slider_max + round(new_range*0.05, 3);
                new_slider_limits = zeros(size(1,2));
                new_slider_limits(1,1) = new_slider_min;
                new_slider_limits(1,2) = new_slider_max;
                E.ui_slider.Limits = new_slider_limits;
            end
            if E.mask_logic ~= E.ui_logic.SelectedObject
                E.mask_logic = E.ui_logic.SelectedObject;
                return_modified = true;
                modified = true;
            end
            if E.mask_cond ~= E.ui_cond.SelectedObject
                E.mask_cond = E.ui_cond.SelectedObject;
                modified = true;
            end
            if E.mask_slider_value ~= E.ui_slider.Value
                E.mask_slider_value = E.ui_slider.Value;
                modified = true;
            end
            if modified
                if E.ui_cond.SelectedObject == E.ui_cond_gt
                    new_mask_array = F.roi & (E.pixvals_smooth > E.mask_slider_value);
                else
                    new_mask_array = F.roi & (E.pixvals_smooth < E.mask_slider_value);
                end
                if isfield(E, 'mask_array') && isequal(new_mask_array, E.mask_array)
                    modified = false;
                else
                    E.mask_array = new_mask_array;
                end
                F.(item_name) = E;
                guidata(fig, F);
            end
        elseif modified
            E.mask_array = F.roi & E.pixvals_smooth;
            F.(item_name) = E;
            guidata(fig, F);
        end
    end
    
    if (modified || redraw) && E.ui_display.Value == 1
        redraw_item(fig, item_num, modified);
    end
    
    modified = modified | return_modified;
end


function recompute_basemask(fig, item_num, modified, redraw)
    if ~exist('modified', 'var') || isempty(modified)
        modified = false;
    end
    if ~exist('redraw', 'var') || isempty(redraw)
        redraw = false;
    end
    
    F = guidata(fig);
    
    if item_num == F.kernel_smooth.item_num
        modified = modified | recompute_basemask_component(fig, F.nans.item_num);
        modified = modified | recompute_basemask_component(fig, F.slope.item_num);
        modified = modified | recompute_basemask_component(fig, F.entropy.item_num);
        modified = modified | recompute_basemask_component(fig, F.elev.item_num);
    else
        modified = modified | recompute_basemask_component(fig, item_num, redraw);
    end
    
    if modified
        F = guidata(fig);
    else
        return;
    end
    
    struct_fields = fieldnames(F);
    for i = 1:numel(struct_fields)
        if (   strcmp(struct_fields{i}, 'slope') ...
            || strcmp(struct_fields{i}, 'entropy') ...
            || strcmp(struct_fields{i}, 'elev'))
            E = F.(struct_fields{i});
            if E.ui_toggle.Value == 1
                if E.ui_logic.SelectedObject == E.ui_logic_or
                    if ~exist('new_basemask_array_or', 'var')
                        new_basemask_array_or = E.mask_array;
                    else
                        new_basemask_array_or = new_basemask_array_or | E.mask_array;
                    end
                else
                    if ~exist('new_basemask_array_and', 'var')
                        new_basemask_array_and = E.mask_array;
                    else
                        new_basemask_array_and = new_basemask_array_and & E.mask_array;
                    end
                end
            end
        end
    end
    
    new_basemask_array = false(size(F.basemask_array));
    if exist('new_basemask_array_or', 'var')
        new_basemask_array = new_basemask_array_or;
    end
    if exist('new_basemask_array_and', 'var')
        if exist('new_basemask_array', 'var')
            new_basemask_array = new_basemask_array | new_basemask_array_and;
        else
            new_basemask_array = new_basemask_array_and;
        end
    end
    if F.nans.ui_toggle.Value == 1
        new_basemask_array = new_basemask_array | F.nans.mask_array;
    end
%     new_basemask_array = new_basemask_array & F.roi;
    
    if isequal(new_basemask_array, F.basemask_array)
        modified = false;
    else
        F.basemask_array = new_basemask_array;
    end
    
    if modified
        fprintf('Call to `recompute_resultmask`\n');
%         recompute_resultmask(fig, item_num, true);
        guidata(fig, F);
%         recompute_resultmask(fig);
    end
end
            



            
            
            
         
    
%     %% Left side %%
%     
%     
%     % Headers
% %     F.left_label_togg = uicontrol('Parent',fig,'Style','text','Position',[260,145,100,23],...
% %                     'String','Kernel Size','BackgroundColor',bgcolor);
%     
%     bgcolor = fig.Color;
%     
%     % Kernel Size
%     kern_min = 1;
%     kern_max = 20;
%     F.kern_val = 5;
%     F.kern_title = uicontrol('Parent',fig,'Style','text','Position',[260,145,100,23],...
%                     'String','Kernel Size','BackgroundColor',bgcolor);
%     F.kern_slide = uicontrol('Parent',fig,'Style','slider','Position',[101,174,419,23],...
%                     'Value',F.kern_val, 'min',kern_min, 'max',kern_max, 'SliderStep',[1/20, 0.1]);
%     F.kern_lmin = uicontrol('Parent',fig,'Style','text','Position',[50,170,43,23],...
%                     'String',num2str(kern_min),'BackgroundColor',bgcolor,...
%                     'HorizontalAlignment','Right');
%     F.kern_lmax = uicontrol('Parent',fig,'Style','text','Position',[528,170,43,23],...
%                     'String',num2str(kern_max),'BackgroundColor',bgcolor,...
%                     'HorizontalAlignment','Left');
%     F.kern_slide.Callback = @(es,ed) callback_kern_slide(fig, es.Value);
% 
% 
%     % Slope Filter
%     F.group_slope = hggroup;
%     F.slope_on = 0;
%     F.kern_val_slope = 0;
%     [dx,dy] = gradient(F.z,F.x,F.y);
%     [~,F.slope_raw] = cart2pol(dx,dy);
%     F.slope_smooth = [];
%     clear dx;
%     clear dy;
% 
%     roi_slope = F.slope_raw;
%     roi_slope(~F.roi) = NaN;
%     slope_min = min(roi_slope(:));
%     slope_max = max(roi_slope(:));
%     F.slope_val = (slope_min+slope_max)/2;
%     clear roi_slope;
% 
%     F.slope_title = uicontrol('Parent',fig,'Style','text','Position',[260,85,100,23],...
%                     'String','Slope Filter','BackgroundColor',bgcolor);
%     F.slope_togg = uicontrol('Parent',fig,'Style','togglebutton','Position',[20,114,23,23],...
%                     'Value',F.slope_on, 'min',0, 'max',1);
%     F.slope_slide = uicontrol('Parent',fig,'Style','slider','Position',[101,114,419,23],...
%                     'Value',F.slope_val, 'min',slope_min, 'max',slope_max);
%     F.slope_lmin = uicontrol('Parent',fig,'Style','text','Position',[50,110,43,23],...
%                     'String',[num2str(round(slope_min*100)),' %'],'BackgroundColor',bgcolor,...
%                     'HorizontalAlignment','Right');
%     F.slope_lmax = uicontrol('Parent',fig,'Style','text','Position',[528,110,43,23],...
%                     'String',[num2str(round(slope_max*100)),' %'],'BackgroundColor',bgcolor,...
%                     'HorizontalAlignment','Left');
%     F.slope_slide.Callback = @(es,ed) callback_slope_slide(fig, es.Value);
%     F.slope_togg.Callback = @(es,ed) callback_slope_togg(fig, es.Value);
% 
% 
%     % Elevation Filter
%     F.group_elev = hggroup;
%     F.elev_on = 0;
%     F.kern_val_elev = 0;
%     F.elev_raw = F.z;
%     F.elev_smooth = [];
% 
%     roi_elev = F.elev_raw;
%     roi_elev(~F.roi) = NaN;
%     elev_min = min(roi_elev(:));
%     elev_max = max(roi_elev(:));
%     F.elev_val = (elev_min+elev_max)/2;
%     clear roi_elev;
% 
%     F.elev_title = uicontrol('Parent',fig,'Style','text','Position',[260,25,100,23],...
%                     'String','Elevation Filter','BackgroundColor',bgcolor);
%     F.elev_togg = uicontrol('Parent',fig,'Style','togglebutton','Position',[20,54,23,23],...
%                     'Value',F.elev_on, 'min',0, 'max',1);
%     F.elev_slide = uicontrol('Parent',fig,'Style','slider','Position',[101,54,419,23],...
%                     'Value',F.elev_val, 'min',elev_min, 'max',elev_max);
%     F.elev_lmin = uicontrol('Parent',fig,'Style','text','Position',[50,50,43,23],...
%                     'String',[num2str(round(elev_min)),' m'],'BackgroundColor',bgcolor,...
%                     'HorizontalAlignment','Right');
%     F.elev_lmax = uicontrol('Parent',fig,'Style','text','Position',[528,50,43,23],...
%                     'String',[num2str(round(elev_max)),' m'],'BackgroundColor',bgcolor,...
%                     'HorizontalAlignment','Left');
%     F.elev_slide.Callback = @(es,ed) callback_elev_slide(fig, es.Value);
%     F.elev_togg.Callback = @(es,ed) callback_elev_togg(fig, es.Value);
% 
% 
%     guidata(fig, F);
%     callback_kern_slide(fig, F.kern_val);
%     F = guidata(fig);
% 
% 
%     prompt = 'apply filters?';
%     while true
%         F.mask_elev = false(size(F.z));
%         F.mask_slope = false(size(F.z));
% 
%         guidata(fig, F);
%         while true
%             s=input([prompt,' (1/0)\n'],'s');
%             if ~strcmpi(s,'1') && ~strcmpi(s,'0')
%                 fprintf('%s not recogized, try again\n',s);
%             else
%                 break
%             end
%         end
%         F = guidata(fig);
%         
%         if strcmpi(s,'1')
%             
%             mask_additions = true(size(F.z));
%             if F.elev_on
%                 mask_additions = (mask_additions & F.mask_elev);
%             end
%             if F.slope_on
%                 mask_additions = (mask_additions & F.mask_slope);
%             end
%             F.mask = F.roi & (F.mask | mask_additions);
% 
%             prompt = 'apply more filters?';
%         else
%             
%             B = bwboundaries(F.mask, 'noholes');
%             polys = [polys; B];
%             
%             break
%         end
%     end
% 
%     delete(F.group_slope);
%     delete(F.group_elev);
% 
%     delete(poly_group);
%     poly_group = hggroup;
%     for k = 1:length(polys)
%        boundary = polys{k};
%        plot(F.x(boundary(:,2)), F.y(boundary(:,1)), 'g','linewidth',2,'Parent',poly_group)
%     end
%     
%     children = get(fig, 'children');
%     for k = 1:numel(children)
%         if strcmp(class(children(k)), 'matlab.ui.control.UIControl')
%             delete(children(k));
%         end
%     end
%     clear F;
% end


function callback_kern_slide(fig, val)
    F = guidata(fig);
    F.kern_val = round(val);
    F.kernel_smooth = ones(F.kern_val)/(F.kern_val.^2);
    guidata(fig, F);
    
    callback_slope_slide(fig, F.slope_val);
    callback_elev_slide(fig, F.elev_val);
end


function callback_slope_togg(fig, val)
    F = guidata(fig);
    F.slope_on = val;
    guidata(fig, F);
    
    callback_slope_slide(fig, F.slope_val)
end


function callback_slope_slide(fig, val)
    F = guidata(fig);
    
    if F.slope_on
        update_mask = false;
        if F.kern_val_slope ~= F.kern_val
            update_mask = true;
            F.kern_val_slope = F.kern_val;
            if F.kern_val > 1
                F.slope_smooth = conv2(F.slope_raw, F.kernel_smooth, 'same');
            else
                F.slope_smooth = F.slope_raw;
            end
        end
        if val ~= F.slope_val
            update_mask = true;
            F.slope_val = val;
        end
        
        if update_mask
%             F.mask_slope = F.roi & ((F.slope_smooth > F.slope_val) | isnan(F.slope_smooth));
            F.mask_slope = F.roi & (F.slope_smooth > F.slope_val);
        end
        
        children = get(F.group_slope, 'children');
        if ~isempty(children)
            delete(children(:));
        end
        B = bwboundaries(F.mask_slope, 'noholes');
        for k = 1:length(B)
           boundary = B{k};
           plot(F.x(boundary(:,2)), F.y(boundary(:,1)), 'magenta','linewidth',2,'Parent',F.group_slope)
        end
    else
        children = get(F.group_slope, 'children');
        if ~isempty(children)
            delete(children(:));
        end
    end
    
    guidata(fig, F);
end


function callback_elev_togg(fig, val)
    F = guidata(fig);
    F.elev_on = val;
    guidata(fig, F);
    
    callback_elev_slide(fig, F.elev_val)
end


function callback_elev_slide(fig, val)
    F = guidata(fig);
    
    if F.elev_on
        update_mask = false;
        if F.kern_val_elev ~= F.kern_val
            update_mask = true;
            F.kern_val_elev = F.kern_val;
            if F.kern_val > 1
                F.elev_smooth = conv2(F.elev_raw, F.kernel_smooth, 'same');
            else
                F.elev_smooth = F.elev_raw;
            end
        end
        if val ~= F.elev_val
            update_mask = true;
            F.elev_val = val;
        end
        
        if update_mask
            F.mask_elev = F.roi & ((F.elev_smooth < F.elev_val) | isnan(F.elev_smooth));
        end
        
        children = get(F.group_elev, 'children');
        if ~isempty(children)
            delete(children(:));
        end
        B = bwboundaries(F.mask_elev, 'noholes');
        for k = 1:length(B)
           boundary = B{k};
           plot(F.x(boundary(:,2)), F.y(boundary(:,1)), 'cyan','linewidth',2,'Parent',F.group_elev)
        end
    else
        children = get(F.group_elev, 'children');
        if ~isempty(children)
            delete(children(:));
        end
    end
    
    guidata(fig, F);
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
