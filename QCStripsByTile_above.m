function [] = QCStripsByTile_above(regionNum,varargin)

global fig_qc
global fig_panel
global filter_item_plot_groups
global filtermode_avail
fig_qc = figure('Name','QCSBT');
fig_panel = [];
filter_item_plot_groups = {};

version_str = sprintf('%s', version);
version_yr = str2num(version_str((end-5):(end-2)));
if version_yr >= 2016
    filtermode_avail = true;
else
    filtermode_avail = false;
    fprintf(2, 'filter mode requires MATLAB release R2016a or later; updating to the latest release is recommended\n');
end

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
qcAll = false;

for i=1:2:length(varargin)
    if strcmp(varargin{i}, 'all')
        qcAll = true;
        continue;
    elseif strcmp(varargin{i}, 'nofilter')
        filtermode_avail = false;
        continue;
    end
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
    for i=1:length(meta.f)
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
    
    qctile(tile,meta,minN,minArea,changePath,qcAll);
    
end
end


function qctile(tiles,meta,minN,minArea,changePath,qcAll)
global fig_qc
global fig_panel
global filter_item_plot_groups
global filtermode_avail

coverageFile_warned = false;
orthoFile_warned = false;
demFile_warned = false;

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


% remove already added entries from metadata
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
    if i>1; for p=1:count; fprintf('\b'); end %delete line before
        count = fprintf('strip %d of %d',i,size(meta.f,1));
    end
    
    coverageFile = strrep(meta.f{i}, 'meta.txt', 'dem_coverage.tif');
    if exist(coverageFile, 'file')
        BW = get_tile_size_coverage(coverageFile, x, y);
        if isempty(BW)
            fprintf(2, ['no overlap between strip and tile; ' ...
                        'make sure *dem_coverage.tif is properly aligned to tile grid\n']);
            continue;
        end
    else
        % locate grid pixels within footprint polygon
        BW = roipoly(x, y, N, meta.x{i}, meta.y{i});
    end
    
    %if mask data exists, apply it
    if meta.qc(i) == 3
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
for p=1:count; fprintf('\b'); end

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
F = [];
while length(meta.f) >= 1
    
    percent_filled = 100*sum(N(:))./numel(N);
    
    if ~qcAll
        if percent_filled == 100; fprintf('100%% filled returning \n',percent_filled); return; end
    end
    
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
    
    if ~qcAll
        
        % remove rendundant files (with already 100% coverage)
        redundantFlag = meta.gridPointN < minN;

        if any(redundantFlag)

            % remove redundant files from lists
            meta = structfun(@(x) ( x(~redundantFlag,:) ), meta, 'UniformOutput', false);

            fprintf('%d redundant files (N < %d) removed\n',sum(redundantFlag),minN);

            if isempty(meta.f); fprintf('all strips removed, returning \n'); return; end
        end
        
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

    hillFile=strrep(meta.f{n},'meta.txt','dem_browse.tif');
    demFile=strrep(meta.f{n},'meta.txt','dem_10m.tif');
    orthoFile=strrep(meta.f{n},'meta.txt','ortho_browse.tif');
    
    fprintf('%s\n',hillFile);
    fprintf('%d new pointsm, gcp sigma=%.2f, mean coreg RMSE=%.2f, max coreg RMSE=%.2f \n',...
        meta.gridPointN(n),meta.sigma(n),meta.avg_rmse(n), meta.max_rmse(n));
    
    I=readGeotiff(hillFile);
    
    
    if isempty(F)
        F = struct();
    else
        F = guidata(fig_qc);
        gui_fields = fieldnames(F);
        for i = 1:length(gui_fields)
            F = rmfield(F, gui_fields{i});
        end
    end
    F.plot_group_hill = hggroup;
    F.plot_group_ortho = hggroup;
    F.fileName_ortho = orthoFile;
    F.ortho_browse_ready = false;
    guidata(fig_qc, F);
    
    % Make square image to show in the figure window.
    [Z_fig, Ni_fig] = get_square_browse_images(fig_qc, I, N, x, y);
    F = guidata(fig_qc);
    
    
    % Plot strip hillshade and tile coverage in figure.
    imagesc(F.X_fig,F.Y_fig,Z_fig,'alphadata',single(Z_fig ~= 0), 'Parent',F.plot_group_hill);
    set(gca,'color','r');
    axis xy  equal tight;
    colormap gray;
    hold on;
    imagesc(F.X_fig,F.Y_fig,Ni_fig,'alphadata',single(Ni_fig).*.25, 'Parent',F.plot_group_hill);
    
    
    % Make radio buttons for switching between hillshade, ortho.
    F.ui_image_select = uibuttongroup('Visible','off',...
                      'Units','pixels',...
                      'Position',[0 0 130 100],...
                      'SelectionChangedFcn',@(source,event) change_browse_image(fig_qc));
    uicontrol(F.ui_image_select,'Style','radiobutton', 'Position',[20 50 100 30],...
                      'String','DEM Hillshade', 'HandleVisibility','off');
    uicontrol(F.ui_image_select,'Style','radiobutton', 'Position',[20 20 100 30],...
                      'String','Ortho Image', 'HandleVisibility','off');
    F.ui_image_select.Visible = 'on';
    guidata(fig_qc, F);
    

    if isfield(tiles,'coastline')
        for i=1:length(tiles.coastline{1})
            if  isempty(tiles.coastline{1}{i}); continue; end
            plot(tiles.coastline{1}{i}(1,:),...
                tiles.coastline{1}{i}(2,:),'b','linewidth',1)
        end
    end
    
    % Plot tile boundary.
    plot([tiles.x0,tiles.x0,tiles.x1,tiles.x1,tiles.x0], [tiles.y0,tiles.y1,tiles.y1,tiles.y0,tiles.y0],'w','linewidth',2)
    
    % Give figure axis a border beyond DEM extent.
    set(gca,'xlim',[min(F.X_fig)-500 max(F.X_fig)+500],'ylim',[min(F.Y_fig)-500 max(F.Y_fig)+500]);
    
    
    %set(gcf,'units','normalized');
    %set(gcf,'position',[0.01,0.01,.35,.9])
    qc=load([fileparts(hillFile),'/qc.mat']);
    
    fileNames = qc.fileNames;
    
   % alter paths in database if set
    if ~isempty(changePath)
        fileNames = strrep(fileNames,'/mnt/pgc',changePath);
        fileNames = strrep(fileNames,'/','\');        
    end
    
    [~,IA]=intersect(fileNames, hillFile);
    
    if isempty(IA) 
        error('this file name not matched in the qc.mat, probably need to upadate it.'); 
    end
    
    if qc.flag(IA) ~= 4 && qc.flag(IA) ~= 0
        
        fprintf('flag previoulsy changed to %d, applying\n',qc.flag(IA))
        
    else
        if ~coverageFile_warned && ~exist(coverageFile, 'file')
            fprintf(2, '** MISSING *_dem_coverage.tif FILE(S); %% FILLED MAY BE INCORRECT UPON QC REVISIT **\n');
            fprintf('(the preceding warning will now be suppressed)\n');
            coverageFile_warned = true;
        end
        if ~orthoFile_warned && ~exist(orthoFile, 'file')
            fprintf(2, '** Missing *_ortho_browse.tif file; ortho image layer is NOT available for this image **\n');
            fprintf('(the preceding warning will now be suppressed)\n');
            orthoFile_warned = true;
        end
        if ~demFile_warned && ~exist(demFile, 'file') && filtermode_avail
            fprintf(2, '** Missing *_dem_10m.tif file; filter mode is NOT available for this image **\n');
            fprintf('(the preceding warning will now be suppressed)\n');
            demFile_warned = true;
        end
        
        while true
            
            try
                
                if qc.flag(IA) == 4
                    fprintf(2, 'qc flag previously set to 4\n');
                end
                    
                flag=input('Enter quality flag: 0=skip,9=back, 1=good, 2=partial, 3=manual edit, 4=poor, 5=unuseable, 6=next tile\n');
                
                
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
        
        
        if flag == 6; clf(fig_qc); return; end
        
        
        if flag == 0;  skipn = skipn+1;  clf(fig_qc); continue; end
        
        if flag == 9;  skipn = skipn-1;  clf(fig_qc); continue; end
        
        
        qc.flag(IA)=flag;
        
        if flag == 1 || flag == 2
            qc.x{IA}=cell(1); qc.y{IA}=cell(1);
            
        end
        
        if flag == 3
            fprintf('entering manual edit mode\n')
            
            total_poly_num = 1;
            while true
                
                [roi_BW,roi_x,roi_y] = roipoly;
                
                if ~isempty(roi_BW)
                    roi_poly = {[roi_x roi_y]};
                    plot_group_roi = plot(roi_x,roi_y,'g','linewidth',2);
                    qc_polys = roi_poly;
                    plot_group_qc = plot_group_roi;

                    if filtermode_avail && exist(demFile, 'file')
                        while true
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
                                I_dem = readGeotiff(demFile); 
                            end
                            if ~exist('I_ortho', 'var')
                                if isfield(F, 'I_ortho')
                                    I_ortho = F.I_ortho;
                                else
                                    if exist(orthoFile, 'file')
                                        I_ortho = readGeotiff(orthoFile);
                                        F.I_ortho = I_ortho;
                                        guidata(fig_qc, F);
                                    else
                                        I_ortho = [];
                                        fprintf(2, '** Missing *_ortho_browse.tif file; falling back to DEM for entropy filter **\n');
                                    end
                                end
                            end
                            
                            % Crop figure-size ROI mask to DEM extent.
                            roi_BW = roi_BW(F.Z_r0:F.Z_r1, F.Z_c0:F.Z_c1);
                            
                            % Enter filter mode.
                            plot_group_roi_black = plot(roi_x,roi_y,'black','linewidth',2);
                            plot_group_filt = hggroup;
                            filt_polys = [];
                            try
                                [filt_polys] = filter_mode(roi_BW, I_dem, I_ortho, plot_group_filt);
                            catch ME
                                if strcmp(ME.identifier, 'MATLAB:guidata:InvalidInput')
                                    fig_panel = [];
                                else
                                    rethrow(ME);
                                end
                            end
                            for i = 1:length(filter_item_plot_groups)
                                delete(filter_item_plot_groups{i});
                            end
                            if ~isempty(filt_polys)
                                qc_polys = filt_polys;
                                plot_group_qc = plot_group_filt;
                                delete(plot_group_roi);
                            else
                                delete(plot_group_filt);
                                delete(plot_group_roi_black);
                            end
                        end
                    end
                    
                    while true
                        s=input('apply mask? (1/0)\n','s');
                        if ~strcmpi(s,'1') && ~strcmpi(s,'0')
                            fprintf('%s not recogized, try again\n',s);
                        else
                            break
                        end
                    end

                    if strcmpi(s,'1')
                        % Write poly(s) to qc file.
                        for poly_num = 1:length(qc_polys)
                            qc.x{IA}{total_poly_num} = qc_polys{poly_num}(:,1);
                            qc.y{IA}{total_poly_num} = qc_polys{poly_num}(:,2);
                            total_poly_num = total_poly_num + 1;
                        end
                    else
                        % Clear drawn polys created with this edit.
                        delete(plot_group_qc);
                    end
                    
                    if exist('plot_group_roi_black', 'var')
                        delete(plot_group_roi_black);
                    end
                end
                
                while true
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
        
        if qc.flag(IA) == 3 && isempty(qc.x{IA}{1})
            fprintf('changing qc flag from 3 to 2\n')
            qc.flag(IA)= 2;
        end
        
        save([fileparts(hillFile),'/qc.mat'],'-struct','qc');
        
    end
    
    clf(fig_qc);
    
    if exist('I_dem', 'var')
        clear I_dem;
    end
    if exist('I_ortho', 'var')
        clear I_ortho;
    end
    
    if qc.flag(IA) > 0 && qc.flag(IA) < 4
        
        M = I.z ~=0;
        
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


function [BW] = get_tile_size_coverage(fileName_matchtag, tile_x, tile_y)
        BW = [];
    
        I = readGeotiff(fileName_matchtag);
        I.z = (I.z ~= 0);

        tile_c0 = 1;
        tile_c1 = numel(tile_x);
        tile_r0 = 1;
        tile_r1 = numel(tile_y);
        strip_c0 = 1;
        strip_c1 = numel(I.x);
        strip_r0 = 1;
        strip_r1 = numel(I.y);

        if I.x(1) < tile_x(1)
            strip_c0 = find(I.x == tile_x(1));
            if isempty(strip_c0)
                return;
            end
        else
            tile_c0 = find(I.x(1) == tile_x);
            if isempty(tile_c0)
                return;
            end
        end
        if I.x(end) > tile_x(end)
            strip_c1 = find(I.x == tile_x(end));
            if isempty(strip_c1)
                return;
            end
        else
            tile_c1 = find(I.x(end) == tile_x);
            if isempty(tile_c1)
                return;
            end
        end
        if I.y(1) > tile_y(1)
            strip_r0 = find(I.y == tile_y(1));
            if isempty(strip_r0)
                return;
            end
        else
            tile_r0 = find(I.y(1) == tile_y);
            if isempty(tile_r0)
                return;
            end
        end
        if I.y(end) < tile_y(end)
            strip_r1 = find(I.y == tile_y(end));
            if isempty(strip_r1)
                return;
            end
        else
            tile_r1 = find(I.y(end) == tile_y);
            if isempty(tile_r1)
                return;
            end
        end

        BW = false(numel(tile_y), numel(tile_x));
        BW(tile_r0:tile_r1, tile_c0:tile_c1) = I.z(strip_r0:strip_r1, strip_c0:strip_c1);
end


function [Z_fig, Ni_fig] = get_square_browse_images(fig, I, N, tile_x, tile_y)
    F = guidata(fig);
    
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
    
    Ni=interp2(tile_x,tile_y(:),N,X_fig,Y_fig(:),'*nearest');
    
    Z_fig = zeros(max(size(Y_fig)), max(size(X_fig)));
    Ni_fig = zeros(max(size(Y_fig)), max(size(X_fig)));
    
    Z_c0 = find(I.x(1) == X_fig);
    Z_c1 = find(I.x(end) == X_fig);
    Z_r0 = find(I.y(1) == Y_fig);
    Z_r1 = find(I.y(end) == Y_fig);
    
    Ni_c0 = find(tile_x(1) == X_fig);
    if isempty(Ni_c0)
        Ni_c0 = 1;
    end
    Ni_c1 = find(tile_x(end) == X_fig);
    if isempty(Ni_c1)
        Ni_c1 = length(X_fig);
    end
    Ni_r0 = find(tile_y(1) == Y_fig);
    if isempty(Ni_r0)
        Ni_r0 = 1;
    end
    Ni_r1 = find(tile_y(end) == Y_fig);
    if isempty(Ni_r1)
        Ni_r1 = length(Y_fig);
    end
    
    Z_fig(Z_r0:Z_r1, Z_c0:Z_c1) = I.z;
    Ni_fig(Ni_r0:Ni_r1, Ni_c0:Ni_c1) = Ni(Ni_r0:Ni_r1, Ni_c0:Ni_c1);
    
    F.X_fig = X_fig;
    F.Y_fig = Y_fig;
    F.Z_r0 = Z_r0;
    F.Z_r1 = Z_r1;
    F.Z_c0 = Z_c0;
    F.Z_c1 = Z_c1;
    
    guidata(fig, F);
end


function change_browse_image(fig)
    F = guidata(fig);
    
    if strcmp(F.ui_image_select.SelectedObject.String, 'DEM Hillshade')
        set(F.plot_group_hill, 'visible','on');
        set(F.plot_group_ortho, 'visible','off');
    elseif strcmp(F.ui_image_select.SelectedObject.String, 'Ortho Image')
        if exist(F.fileName_ortho, 'file')
            if ~F.ortho_browse_ready
                if ~isfield(F, 'I_ortho')
                    F.I_ortho = readGeotiff(F.fileName_ortho);
                end
                O_fig = zeros(max(size(F.Y_fig)), max(size(F.X_fig)));
                O_fig(F.Z_r0:F.Z_r1, F.Z_c0:F.Z_c1) = imadjust(F.I_ortho.z, stretchlim(F.I_ortho.z), []);
                imagesc(F.X_fig,F.Y_fig,O_fig, 'Parent',F.plot_group_ortho);
                uistack(F.plot_group_ortho, 'bottom');
                F.ortho_browse_ready = true;
            end
            set(F.plot_group_hill, 'visible','off');
            set(F.plot_group_ortho, 'visible','on');
        end
    end
    
    guidata(fig, F);
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
            if isfield(F, 'ortho')
                z_uint16 = uint16(F.ortho);
            else
                z_uint16 = cast(F.z, 'uint16');
            end
            z_subtraction =  movmax(z_uint16, entropy_kernel_size) - movmin(z_uint16, entropy_kernel_size);
            A = entropyfilt(z_subtraction, true(entropy_kernel_size));
        case F.elev.item_num
            A = F.z;
        otherwise
            error('`get_pixvals_raw` item_num=%d not recognized\n', item_num);
    end
end


function [E, item_name] = get_gui_element_by_item_num(F, item_num)
    gui_fields = fieldnames(F);
    for i = 1:length(gui_fields)
        E = F.(gui_fields{i});
        if isstruct(E)
            element_fields = fieldnames(E);
            for j = 1:length(element_fields)
                if strcmp(element_fields{j}, 'item_num')
                    if E.item_num == item_num
                        item_name = gui_fields{i};
                        return
                    end
                end
            end
        end
    end
    error('`get_gui_element_by_item_num` unsuccessfull` for item_num=%d', item_num);
end


function [polys] = filter_mode(roi_BW, I_dem, I_ortho, plot_group_filt)
    global fig_panel
    global filter_item_plot_groups

    polys = {};
    plot_group_filt_single = hggroup;
    
    % Create Filter Controls UI in a new window.
    if isempty(fig_panel)
        create_filter_panel();
        initial_recompute = false;
    else
        reset_filter_data();
        initial_recompute = true;
    end
    F = guidata(fig_panel);
    gui_fields = fieldnames(F);
    for i = 1:length(gui_fields)
        E = F.(gui_fields{i});
        if isstruct(E) && isfield(E, 'plot_group')
            filter_item_plot_groups(end+1) = {E.plot_group};
        end
    end
    
    % Crop ROI and DEM to ROI extent.
    F = guidata(fig_panel);
    [roi_BW,rcrop,ccrop] = cropzeros(roi_BW);
    F.z = I_dem.z(rcrop(1):rcrop(2), ccrop(1):ccrop(2));
    F.x = I_dem.x(ccrop(1):ccrop(2));
    F.y = I_dem.y(rcrop(1):rcrop(2));
    clear I_dem;
    F.z(F.z == -9999) = NaN;
    if ~isempty(I_ortho) 
        F.ortho = I_ortho.z(rcrop(1):rcrop(2), ccrop(1):ccrop(2));
        clear I_ortho;
    end
    F.roi = roi_BW;
    F.basemask_array = false(size(F.z));
    guidata(fig_panel, F);
    
    if initial_recompute
        recompute_basemask(fig_panel, F.kernel_smooth.item_num, true);
    end
    
    prompt = 'apply filters?';
    quit = false;
    while ~quit
        while true
            s=input([prompt,' (1/0)\n'],'s');
            if ~strcmpi(s,'1') && ~strcmpi(s,'0')
                fprintf('%s not recogized, try again\n',s);
            else
                break
            end
        end
        
        if strcmpi(s,'1')
            F = guidata(fig_panel);
            if F.result.ui_display.Value == 0
                fprintf('Please set RESULT MASK display to SHOW before proceeding.\n');
                continue;
            else
                mask_polys = F.result.mask_polys;
                if ~isempty(mask_polys)
                    new_polys = cellfun(@(p) vertcat(F.x(p(:,2)), F.y(p(:,1))).', mask_polys, 'UniformOutput',false);
                    for k = 1:length(new_polys)
                        p = new_polys{k};
                        plot(p(:,1), p(:,2), 'white','linewidth',2, 'Parent',plot_group_filt_single);
                    end
                    polys = [polys; new_polys];
                end
            end
            prompt = 'apply more filters?';
        else
            quit = true;
        end
    end
    
    delete(plot_group_filt_single);
    for k = 1:length(polys)
        p = polys{k};
        plot(p(:,1), p(:,2), 'green','linewidth',2, 'Parent',plot_group_filt);
    end
    
    reset_filter_data();
end


function reset_filter_data()
    global fig_panel
    
    F = guidata(fig_panel);
    gui_fields = fieldnames(F);
    for i = 1:length(gui_fields)
        E = F.(gui_fields{i});
        if isstruct(E)
            gui_item_fields = fieldnames(E);
            save_fieldnames = {};
            j = 1;
            while j <= length(gui_item_fields)
                fname = gui_item_fields{j};
                if ~isempty(save_fieldnames)
                    save = false;
                    for k = 1:length(save_fieldnames)
                        if strcmp(fname, save_fieldnames{k})
                            save = true;
                            break;
                        end
                    end
                    if ~save
                        E = rmfield(E, fname);
                    end
                elseif strcmp(fname, 'initial_fieldnames')
                    save_fieldnames = E.initial_fieldnames;
                    j = 0;
                    if isfield(E, 'plot_group')
                        delete(E.plot_group);
                        E.plot_group = hggroup;
                    end
                end
                j = j + 1;
            end
            if ~isempty(save_fieldnames)
                F.(gui_fields{i}) = E;
            end
        end
    end
    guidata(fig_panel, F);
end


function toggle_item_on(fig, item_num, btn, isKernel)
    if ~exist('isKernel', 'var') || isempty(isKernel)
            isKernel = false;
    end
    
    F = guidata(fig);
    [E,item_name] = get_gui_element_by_item_num(F, item_num);
    
    if item_num == F.cluster_concav.item_num && btn.Value == 0 && F.cluster_dist.ui_toggle.Value == 1
        btn.Value = 1;
        return;
    end
    
    if btn.Value == 0
        if item_num == F.result.item_num
            btn.Text = 'NORMAL';
        else
            btn.Text = 'OFF';
        end
    else
        if item_num == F.result.item_num
            btn.Text = 'INVERTED';
        else
            btn.Text = 'ON';
        end
    end
    
    if item_num == F.cluster_dist.item_num && btn.Value == 1 && F.cluster_concav.ui_toggle.Value == 0
        F.cluster_concav.ui_toggle.Value = 1;
        F.cluster_concav.ui_toggle.Text = 'ON';
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
    if btn.Value == 1 && item_num == F.result.item_num
        recompute_resultmask(fig, item_num, false, true);
    end
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
    [E,item_name] = get_gui_element_by_item_num(F, item_num);
    
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
    [E,item_name] = get_gui_element_by_item_num(F, item_num);
    
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
    [E,item_name] = get_gui_element_by_item_num(F, item_num);
    
    if ~isfield(E, 'plot_group')
        return;
    end
    
    % Clear drawn polygons for this item.
    prev_plot_polys = get(E.plot_group, 'children');
    if ~isempty(prev_plot_polys)
        delete(prev_plot_polys(:));
    end
    
    if (E.item_num == F.result.item_num || E.ui_toggle.Value == 1) && E.ui_display.Value == 1
        % Draw polygons for this item.
        if modified || ~isfield(E, 'mask_polys')
            E.mask_polys = bwboundaries(E.mask_array, 'noholes');
        end
        for k = 1:length(E.mask_polys)
           poly = E.mask_polys{k};
           plot(F.x(poly(:,2)), F.y(poly(:,1)), E.plot_color,'linewidth',2, 'Parent',E.plot_group);
        end
    end
    
    F.(item_name) = E;
    guidata(fig, F);
end


function [modified] = cluster_filter_size(fig, basemask_modified, slider_modified)
    F = guidata(fig);
    E = F.cluster_size;
    
    modified = basemask_modified | slider_modified;
    
    latest_postfilt_array = F.basemask_array;
    
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
    if ~slider_modified && ~basemask_modified && ~isequal(latest_postfilt_array, E.prefilt_array)
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
            new_slider_limits = zeros(1, 2);
            new_slider_limits(1,1) = new_slider_min;
            new_slider_limits(1,2) = new_slider_max + 0.001*(new_slider_max-new_slider_min);
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
            new_postfilt_array = bwareaopen(E.prefilt_array, round(E.ui_slider.Value));
            if E.ui_cond.SelectedObject == E.ui_cond_lt
                new_postfilt_array = xor(E.prefilt_array, new_postfilt_array);
            end
            E.postfilt_array = new_postfilt_array;
        end
        F.cluster_size = E;
        guidata(fig, F);
    end
end


function [modified] = cluster_dilate(fig, prefilt_modified, slider_modified)
    F = guidata(fig);
    E = F.kernel_dilate;
    
    modified = prefilt_modified | slider_modified;
    
    if F.cluster_size.ui_toggle.Value == 1
        latest_postfilt_array = F.cluster_size.postfilt_array;
    else
        latest_postfilt_array = F.basemask_array;
    end
    
    if ~isfield(E, 'prefilt_array')
        prefilt_modified = true;
        modified = true;
    end
    if E.mask_slider_value ~= E.ui_slider.Value
        E.mask_slider_value = E.ui_slider.Value;
        modified = true;
    end
    if ~slider_modified && ~prefilt_modified && ~isequal(E.prefilt_array, latest_postfilt_array)
        prefilt_modified = true;
        modified = true;
    end
    if modified
        if prefilt_modified
            E.prefilt_array = latest_postfilt_array;
        end
        if E.ui_slider.Value > 1
            E.postfilt_array = F.roi & imdilate(E.prefilt_array, true(E.ui_slider.Value));
        else
            E.postfilt_array = E.prefilt_array;
        end
        F.kernel_dilate = E;
        guidata(fig, F);
    end
end


function [int_coords] = image_coords_float_to_int(float_coords, img_nrows, img_ncols)
    int_coords = round(float_coords);
    int_coords(int_coords == 0) = 1;
    int_coords((int_coords(:,1) > img_ncols), 1) = img_ncols;
    int_coords((int_coords(:,2) > img_nrows), 2) = img_nrows;
end


function [modified] = cluster_filter_dist(fig, basemask_modified, slider_modified, dilated)
    F = guidata(fig);
    E = F.cluster_dist;
    
    modified = basemask_modified | slider_modified | dilated;
    sizefilt_modified = false;
    
    if F.kernel_dilate.ui_toggle.Value == 1
        latest_postfilt_array = F.kernel_dilate.postfilt_array;
    elseif F.cluster_size.ui_toggle.Value == 1
        latest_postfilt_array = F.cluster_size.postfilt_array;
    else
        latest_postfilt_array = F.basemask_array;
    end
    
    if F.cluster_size.ui_toggle.Value == 1
        latest_postfilt_CC_ind = F.cluster_size.postfilt_CC_ind;
    else
        latest_postfilt_CC_ind = [1:F.basemask_CC.NumObjects];
    end
    if ~isfield(E, 'prefilt_CC_ind') || ~isequal(E.prefilt_CC_ind, latest_postfilt_CC_ind)
        E.prefilt_CC_ind = latest_postfilt_CC_ind;
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
    if modified
        if basemask_modified
            rp = regionprops(F.basemask_CC, 'Centroid');
            num_clusters = F.basemask_CC.NumObjects;
            [num_rows, num_cols] = size(latest_postfilt_array);
            cent_coords = image_coords_float_to_int(vertcat(rp.Centroid), num_rows, num_cols);
            F.basemask_RP.cent_coords = cent_coords;
            cent_coords_col = cent_coords(:,1);
            cent_coords_row = cent_coords(:,2);
            F.basemask_RP.dist_array = sqrt( ...
                power(repmat(cent_coords_row, [1,num_clusters]) - repmat(cent_coords_row.', [num_clusters,1]), 2) ...
              + power(repmat(cent_coords_col, [1,num_clusters]) - repmat(cent_coords_col.', [num_clusters,1]), 2));
            sizefilt_modified = true;
        end
        
        if sizefilt_modified
            E.dist_array = F.basemask_RP.dist_array(E.prefilt_CC_ind, E.prefilt_CC_ind);
            new_slider_min = min(F.basemask_RP.dist_array(:));
            new_slider_max = max(F.basemask_RP.dist_array(:));
            if new_slider_min == new_slider_max
                new_slider_max = new_slider_min + 1;
            end
            new_slider_limits = zeros(1, 2);
            new_slider_limits(1,1) = new_slider_min;
            new_slider_limits(1,2) = new_slider_max;
            E.ui_slider.Limits = new_slider_limits;
        end
        
        if ~dilated || ~isfield(E, 'postfilt_CC_ind')
            num_clusters = length(E.prefilt_CC_ind);
            dist_array = E.dist_array;
            dist_filter_value = E.mask_slider_value;
            prefilt_CC_ind = E.prefilt_CC_ind;
            new_postfilt_CC_ind = [];
            for i=1:num_clusters
                for j=1:(i-1)
                    if dist_array(i,j) < dist_filter_value
                        new_postfilt_CC_ind = [new_postfilt_CC_ind; [prefilt_CC_ind(i) prefilt_CC_ind(j)]];
                    end
                end
            end
            if isfield(E, 'postfilt_CC_ind') && isequal(new_postfilt_CC_ind, E.postfilt_CC_ind)
                if slider_modified && ~basemask_modified
                    modified = false;
                elseif isfield(E, 'prefilt_array') && isequal(latest_postfilt_array, E.prefilt_array)
                    modified = false;
                end
            else
                E.postfilt_CC_ind = new_postfilt_CC_ind;
            end
        end
        
        if modified
            if ~isfield(E, 'prefilt_array') || ~slider_modified
                E.prefilt_array = latest_postfilt_array;
            end
            new_postfilt_array = E.prefilt_array;
            array_size = size(new_postfilt_array);
            postfilt_CC_ind = E.postfilt_CC_ind;
            cent_coords = F.basemask_RP.cent_coords;
            for k=1:size(postfilt_CC_ind,1)
                cluster_pair_ind = postfilt_CC_ind(k,:);
                i = cluster_pair_ind(1);
                j = cluster_pair_ind(2);
                x = [cent_coords(i,1) cent_coords(j,1)];
                y = [cent_coords(i,2) cent_coords(j,2)];
                nPoints = max(abs(diff(x)), abs(diff(y))) + 1;
                rIndices = round(linspace(y(1), y(2), nPoints));
                cIndices = round(linspace(x(1), x(2), nPoints));
                indices = sub2ind(array_size, rIndices, cIndices);
                new_postfilt_array(indices) = true;
            end
            E.postfilt_array = new_postfilt_array;
        end
        
        F.cluster_dist = E;
        guidata(fig, F);
    end
end 


function [modified] = cluster_adjust_concavity(fig, prefilt_modified, slider_modified)
    F = guidata(fig);
    E = F.cluster_concav;
    
    modified = prefilt_modified | slider_modified;
    
    if F.cluster_dist.ui_toggle.Value == 1
        latest_item_name = 'cluster_dist';
        latest_postfilt_array = F.cluster_dist.postfilt_array;
    elseif F.kernel_dilate.ui_toggle.Value == 1
        latest_item_name = 'kernel_dilate';
        latest_postfilt_array = F.kernel_dilate.postfilt_array;
    elseif F.cluster_size.ui_toggle.Value == 1
        latest_item_name = 'cluster_size';
        latest_postfilt_array = F.cluster_size.postfilt_array;
    else
        latest_item_name = 'basemask_RP';
        latest_postfilt_array = F.basemask_array;
    end
    
    if ~isfield(E, 'prefilt_array')
        prefilt_modified = true;
        modified = true;
    end
    if E.mask_slider_value ~= E.ui_slider.Value
        E.mask_slider_value = E.ui_slider.Value;
        modified = true;
    end
    if ~slider_modified && ~prefilt_modified && ~isequal(E.prefilt_array, latest_postfilt_array)
        prefilt_modified = true;
        modified = true;
    end
    if modified
        if prefilt_modified
            E.prefilt_array = latest_postfilt_array;
        end
        new_postfilt_polys = [];
%         if E.ui_slider.Value == 0
%             [num_rows, num_cols] = size(latest_postfilt_array);
%             if F.cluster_dist.ui_toggle.Value == 1
%                 if ~isfield(F.cluster_dist, 'postfilt_ConvexHull') || prefilt_modified || slider_modified
%                     cc = bwconncomp(F.cluster_dist.postfilt_array);
%                     rp = regionprops(cc, 'ConvexHull');
%                     F.cluster_dist.postfilt_ConvexHull = cellfun(@(y) image_coords_float_to_int(y, num_cols, num_rows),...
%                         arrayfun(@(x) circshift(x.ConvexHull,1,2), rp, 'UniformOutput',false),...
%                         'UniformOutput',false);
%                 end
%                 new_postfilt_polys = F.cluster_dist.postfilt_ConvexHull;
%             elseif F.kernel_dilate.ui_toggle.Value == 1
%                 if ~isfield(F.kernel_dilate, 'postfilt_ConvexHull') || prefilt_modified || slider_modified
%                     cc = bwconncomp(F.kernel_dilate.postfilt_array);
%                     rp = regionprops(cc, 'ConvexHull');
%                     F.kernel_dilate.postfilt_ConvexHull = cellfun(@(y) image_coords_float_to_int(y, num_cols, num_rows),...
%                         arrayfun(@(x) circshift(x.ConvexHull,1,2), rp, 'UniformOutput',false),...
%                         'UniformOutput',false);
%                 end
%                 new_postfilt_polys = F.kernel_dilate.postfilt_ConvexHull;
%             else
%                 if ~isfield(F.basemask_RP, 'ConvexHull') || prefilt_modified || slider_modified
%                     rp = regionprops(F.basemask_CC, 'ConvexHull');
%                     F.basemask_RP.ConvexHull = cellfun(@(y) image_coords_float_to_int(y, num_cols, num_rows),...
%                         arrayfun(@(x) circshift(x.ConvexHull,1,2), rp, 'UniformOutput',false),...
%                         'UniformOutput',false);
%                 end
%                 if F.cluster_size.ui_toggle.Value == 1
%                     new_postfilt_polys = F.basemask_RP.ConvexHull(F.cluster_size.postfilt_CC_ind);
%                 else
%                     new_postfilt_polys = F.basemask_RP.ConvexHull;
%                 end
%             end
%         else
            latest_item = F.(latest_item_name);
            if ~isfield(latest_item, 'postfilt_bwboundaries') || prefilt_modified
                latest_item.postfilt_bwboundaries = bwboundaries(latest_postfilt_array, 'noholes');
                F.(latest_item_name) = latest_item;
            end
            B = latest_item.postfilt_bwboundaries;
            num_polys = length(B);
            K = cellfun(@(x) boundary(x(:,2),x(:,1),E.ui_slider.Value), B, 'UniformOutput',false);
            new_postfilt_polys = cell(num_polys, 1);
            for i=1:num_polys
                k = K{i};
                new_postfilt_polys{i} = B{i}(k,:);
            end
%         end
        E.postfilt_polys = new_postfilt_polys;
        F.cluster_concav = E;
        guidata(fig, F);
    end
end


function [polys] = get_resultmask_polys(fig)
    F = guidata(fig);
    
    polys = {};
    invert_result = (F.result.ui_toggle.Value == 1);
    
    if ~any(F.basemask_array(:))
        if invert_result
            polys = bwboundaries(F.roi);
        end
        return;
    end
    
    if F.cluster_concav.ui_toggle.Value == 1
        polys = F.cluster_concav.postfilt_polys;
        if invert_result
            latest_postfilt_array = false(size(F.roi));
            for k = 1:length(polys)
                p = polys{k};
                M = poly2mask(p(:,2),p(:,1), size(F.roi,1),size(F.roi,2));
                latest_postfilt_array = latest_postfilt_array | M;
            end
            polys = bwboundaries(~latest_postfilt_array & F.roi, 'noholes');
        end
    else
        if F.kernel_dilate.ui_toggle.Value == 1
            latest_item_name = 'kernel_dilate';
            latest_postfilt_array = F.kernel_dilate.postfilt_array;
        elseif F.cluster_size.ui_toggle.Value == 1
            latest_item_name = 'cluster_size';
            latest_postfilt_array = F.cluster_size.postfilt_array;
        else
            latest_item_name = 'basemask_RP';
            latest_postfilt_array = F.basemask_array;
        end
        
        latest_item = F.(latest_item_name);
        latest_item.postfilt_bwboundaries = bwboundaries(latest_postfilt_array, 'noholes');
        F.(latest_item_name) = latest_item;
        if invert_result
            polys = bwboundaries(~latest_postfilt_array & F.roi, 'noholes');
        else
            polys = latest_item.postfilt_bwboundaries;
        end
        
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
            || F.cluster_dist.ui_toggle.Value == 1 ...
            || F.cluster_concav.ui_toggle.Value == 1) ...
        && (~isfield(F, 'basemask_CC') || item_num < F.cluster_size.item_num))
        F.basemask_CC = bwconncomp(F.basemask_array);
        F.basemask_RP = struct();
        guidata(fig, F);
    elseif item_num < F.cluster_size.item_num
        if isfield(F, 'basemask_CC')
            F = rmfield(F, 'basemask_CC');
        end
        F.basemask_RP = struct();
        guidata(fig, F);
    end
    
    if item_num == F.kernel_dilate.item_num
        prev_cluster_dilate_size = F.kernel_dilate.cluster_dilate_size;
        if F.kernel_dilate.ui_toggle.Value == 1
            F.kernel_dilate.cluster_dilate_size = F.kernel_dilate.ui_slider.Value;
        else
            F.kernel_dilate.cluster_dilate_size = 1;
        end
        guidata(fig, F);
        if prev_cluster_dilate_size == F.kernel_dilate.cluster_dilate_size
            return;
        end
    end
    
    clear_result = false;
    if any(F.basemask_array(:))
        if F.cluster_size.ui_toggle.Value == 1 && (item_num == F.cluster_size.item_num || (item_num < F.cluster_size.item_num && modified))
            modified = cluster_filter_size(fig, modified, (cluster_slider && item_num == F.cluster_size.item_num));
        end
        if item_num == F.cluster_size.item_num && redraw
            modified = true;
        end
        if F.kernel_dilate.ui_toggle.Value == 1 && (item_num == F.kernel_dilate.item_num || (item_num < F.kernel_dilate.item_num && modified))
            modified = cluster_dilate(fig, modified, (cluster_slider && item_num == F.kernel_dilate.item_num));
        end
        if item_num == F.kernel_dilate.item_num && redraw
            modified = true;
        end
        if F.cluster_dist.ui_toggle.Value == 1 && (item_num == F.cluster_dist.item_num || (item_num < F.cluster_dist.item_num && modified))
            modified = cluster_filter_dist(fig, (modified && item_num < F.cluster_size.item_num), (cluster_slider && item_num == F.cluster_dist.item_num), item_num == F.kernel_dilate.item_num);
        end
        if item_num == F.cluster_dist.item_num && redraw
            modified = true;
        end
        if F.cluster_concav.ui_toggle.Value == 1 && (item_num == F.cluster_concav.item_num || (item_num < F.cluster_concav.item_num && modified))
            modified = cluster_adjust_concavity(fig, modified, (cluster_slider && item_num == F.cluster_concav.item_num));
        end
        if item_num == F.cluster_concav.item_num && redraw
            modified = true;
        end
    else
        clear_result = true;
    end
    
    if (~isfield(F.result, 'mask_polys') || modified || redraw) && F.result.ui_display.Value == 1
        if ~clear_result
            mask_polys = get_resultmask_polys(fig);
        else
            mask_polys = {};
        end
        F = guidata(fig);
        
        F.result.mask_polys = mask_polys;
        
        guidata(fig, F);
        redraw_item(fig, F.result.item_num, false);
    end
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
    for i = 1:length(struct_fields)
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
        guidata(fig, F);
%         fprintf('Call to `recompute_resultmask`\n');
        recompute_resultmask(fig, item_num, true, true);
    end
end


function [modified] = recompute_basemask_component(fig, item_num, redraw)
    if ~exist('redraw', 'var') || isempty(redraw)
        redraw = false;
    end

    F = guidata(fig);
    [E, item_name] = get_gui_element_by_item_num(F, item_num);
    
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
                new_slider_limits = zeros(1, 2);
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


function create_filter_panel()
    global fig_panel
    
    fig = uifigure('Name','Filter Controls', 'Position',[0 0 960 410]);
    F = guihandles(fig);
    
    % NANS %
    E = struct();
    E.item_num = 1;
    E.item_group = 1;
    E.plot_group = hggroup;
    E.plot_color = 'magenta';
    E.mask_kernel_size = 1;
    E.ui_title = uilabel(fig, 'Position',[155 362 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','Mask NaNs');
    E.ui_toggle = uibutton(fig, 'state', 'Position',[115 355 35 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_on(fig, E.item_num, btn, false),...
        'Value',false, 'Text','OFF');  % DEFAULT
    E.ui_display = uibutton(fig, 'state', 'Position',[430 355 50 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_display(fig, E.item_num, btn),...
        'Value',true, 'Text','SHOW');  % DEFAULT
    E.initial_fieldnames = [];
    E.initial_fieldnames = fieldnames(E);
    F.nans = E;
    
    % SMOOTH KERNEL SIZE %
    E = struct();
    E.item_num = 2;
    E.item_group = 1;
    max_kernel_size = 30;  % SETTING
    E.ui_title = uilabel(fig, 'Position',[155 300 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','Smooth Kernel Size (sq. side length)');
    E.ui_toggle = uibutton(fig, 'state', 'Position',[115 275 35 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_on(fig, E.item_num, btn, true),...
        'Value',false, 'Text','OFF');  % DEFAULT
    E.ui_slider_dec = uibutton(fig, 'push', 'Position',[160 280 20 25], 'Text','-',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, -1, true, true));
    E.ui_slider = uislider(fig, 'Position',[190 290 200 3],...
        'ValueChangedFcn',@(sld,event) slider_drag_integer(fig, E.item_num),...
        'MajorTicksMode','manual', 'MajorTicks',[0:5:max_kernel_size],...
        'MinorTicksMode','manual', 'MinorTicks',[1:1:max_kernel_size],...
        'Limits',[1 max_kernel_size], 'Value',5);  % DEFAULT
    E.ui_slider_inc = uibutton(fig, 'push', 'Position',[400 280 20 25], 'Text','+',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, 1, true, true));
    E.initial_fieldnames = [];
    E.initial_fieldnames = fieldnames(E);
    F.kernel_smooth = E;
    
    % SLOPE FILTER %
    E = struct();
    E.item_num = 3;
    E.item_group = 1;
    E.plot_group = hggroup;
    E.plot_color = 'cyan';
    E.mask_kernel_size = 1;
    E.ui_title = uilabel(fig, 'Position',[155 220 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','Slope Filter (% grade)');
    E.ui_logic = uibuttongroup(fig, 'Position',[20 192 50 40],...
        'SelectionChangedFcn',@(bg,event) recompute_basemask(fig, E.item_num));
    E.ui_logic_or = uiradiobutton(E.ui_logic, 'Position',[3 7 45 30], 'Text','OR');
    E.ui_logic_and = uiradiobutton(E.ui_logic, 'Position',[3 -11 45 30], 'Text','AND');
    E.ui_logic.SelectedObject = E.ui_logic_or;  % DEFAULT
    E.mask_logic = E.ui_logic.SelectedObject;
    E.ui_cond = uibuttongroup(fig, 'Position',[70 192 35 40],...
        'SelectionChangedFcn',@(bg,event) recompute_basemask(fig, E.item_num));
    E.ui_cond_gt = uiradiobutton(E.ui_cond, 'Position',[5 7 30 30], 'Text','>');
    E.ui_cond_lt = uiradiobutton(E.ui_cond, 'Position',[5 -11 30 30], 'Text','<');
    E.ui_cond.SelectedObject = E.ui_cond_gt;  % DEFAULT
    E.mask_cond = E.ui_cond.SelectedObject;
    E.ui_toggle = uibutton(fig, 'state', 'Position',[115 195 35 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_on(fig, E.item_num, btn, false),...
        'Value',false, 'Text','OFF');  % DEFAULT
    E.ui_slider_dec = uibutton(fig, 'push', 'Position',[160 200 20 25], 'Text','-',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, -1, false));
    E.ui_slider = uislider(fig, 'Position',[190 210 200 3],...
        'ValueChangedFcn',@(sld,event) recompute_basemask(fig, E.item_num),...
        'Limits',[0 100], 'Value',0);  % DEFAULT
    E.ui_slider_inc = uibutton(fig, 'push', 'Position',[400 200 20 25], 'Text','+',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, 1, false));
    E.mask_slider_value = E.ui_slider.Value;
    E.ui_display = uibutton(fig, 'state', 'Position',[430 195 50 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_display(fig, E.item_num, btn),...
        'Value',true, 'Text','SHOW');  % DEFAULT
    E.initial_fieldnames = [];
    E.initial_fieldnames = fieldnames(E);
    F.slope = E;
    
    % ENTROPY FILTER %
    E = struct();
    E.item_num = 4;
    E.item_group = 1;
    E.plot_group = hggroup;
    E.plot_color = 'blue';
    E.mask_kernel_size = 1;
    E.ui_title = uilabel(fig, 'Position',[155 140 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','Entropy Filter');
    E.ui_logic = uibuttongroup(fig, 'Position',[20 112 50 40],...
        'SelectionChangedFcn',@(bg,event) recompute_basemask(fig, E.item_num));
    E.ui_logic_or = uiradiobutton(E.ui_logic, 'Position',[3 7 45 30], 'Text','OR');
    E.ui_logic_and = uiradiobutton(E.ui_logic, 'Position',[3 -11 45 30], 'Text','AND');
    E.ui_logic.SelectedObject = E.ui_logic_or;  % DEFAULT
    E.mask_logic = E.ui_logic.SelectedObject;
    E.ui_cond = uibuttongroup(fig, 'Position',[70 112 35 40],...
        'SelectionChangedFcn',@(bg,event) recompute_basemask(fig, E.item_num));
    E.ui_cond_gt = uiradiobutton(E.ui_cond, 'Position',[5 7 30 30], 'Text','>');
    E.ui_cond_lt = uiradiobutton(E.ui_cond, 'Position',[5 -11 30 30], 'Text','<');
    E.ui_cond.SelectedObject = E.ui_cond_gt;  % DEFAULT
    E.mask_cond = E.ui_cond.SelectedObject;
    E.ui_toggle = uibutton(fig, 'state', 'Position',[115 115 35 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_on(fig, E.item_num, btn, false),...
        'Value',false, 'Text','OFF');  % DEFAULT
    E.ui_slider_dec = uibutton(fig, 'push', 'Position',[160 120 20 25], 'Text','-',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, -1, false));
    E.ui_slider = uislider(fig, 'Position',[190 130 200 3],...
        'ValueChangedFcn',@(sld,event) recompute_basemask(fig, E.item_num),...
        'Limits',[0 10], 'Value',0);  % DEFAULT
    E.ui_slider_inc = uibutton(fig, 'push', 'Position',[400 120 20 25], 'Text','+',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, 1, false));
    E.mask_slider_value = E.ui_slider.Value;
    E.ui_display = uibutton(fig, 'state', 'Position',[430 115 50 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_display(fig, E.item_num, btn),...
        'Value',true, 'Text','SHOW');  % DEFAULT
    E.initial_fieldnames = [];
    E.initial_fieldnames = fieldnames(E);
    F.entropy = E;
    
    % ELEVATION FILTER %
    E = struct();
    E.item_num = 5;
    E.item_group = 1;
    E.plot_group = hggroup;
    E.plot_color = 'yellow';
    E.mask_kernel_size = 1;
    E.ui_title = uilabel(fig, 'Position',[155 60 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','Elevation Filter (m)');
    E.ui_logic = uibuttongroup(fig, 'Position',[20 32 50 40],...
        'SelectionChangedFcn',@(bg,event) recompute_basemask(fig, E.item_num));
    E.ui_logic_or = uiradiobutton(E.ui_logic, 'Position',[3 7 45 30], 'Text','OR');
    E.ui_logic_and = uiradiobutton(E.ui_logic, 'Position',[3 -11 45 30], 'Text','AND');
    E.ui_logic.SelectedObject = E.ui_logic_or;  % DEFAULT
    E.mask_logic = E.ui_logic.SelectedObject;
    E.ui_cond = uibuttongroup(fig, 'Position',[70 32 35 40],...
        'SelectionChangedFcn',@(bg,event) recompute_basemask(fig, E.item_num));
    E.ui_cond_gt = uiradiobutton(E.ui_cond, 'Position',[5 7 30 30], 'Text','>');
    E.ui_cond_lt = uiradiobutton(E.ui_cond, 'Position',[5 -11 30 30], 'Text','<');
    E.ui_cond.SelectedObject = E.ui_cond_gt;  % DEFAULT
    E.mask_cond = E.ui_cond.SelectedObject;
    E.ui_toggle = uibutton(fig, 'state', 'Position',[115 35 35 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_on(fig, E.item_num, btn, false),...
        'Value',false, 'Text','OFF');  % DEFAULT
    E.ui_slider_dec = uibutton(fig, 'push', 'Position',[160 40 20 25], 'Text','-',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, -1, false));
    E.ui_slider = uislider(fig, 'Position',[190 50 200 3],...
        'ValueChangedFcn',@(sld,event) recompute_basemask(fig, E.item_num),...
        'Limits',[0 100], 'Value',0);  % DEFAULT
    E.ui_slider_inc = uibutton(fig, 'push', 'Position',[400 40 20 25], 'Text','+',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, 1, false));
    E.mask_slider_value = E.ui_slider.Value;
    E.ui_display = uibutton(fig, 'state', 'Position',[430 35 50 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_display(fig, E.item_num, btn),...
        'Value',true, 'Text','SHOW');  % DEFAULT
    E.initial_fieldnames = [];
    E.initial_fieldnames = fieldnames(E);
    F.elev = E;
    
    
    % RESULT MASK %
    E = struct();
    E.item_num = 10;
    E.item_group = 2;
    E.plot_group = hggroup;
    E.plot_color = 'green';
    E.ui_title = uilabel(fig, 'Position',[615 362 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','RESULT MASK');
    E.ui_toggle = uibutton(fig, 'state', 'Position',[530 355 80 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_on(fig, E.item_num, btn, false),...
        'Value',false, 'Text','NORMAL');  % DEFAULT
    E.ui_display = uibutton(fig, 'state', 'Position',[890 355 50 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_display(fig, E.item_num, btn),...
        'Value',true, 'Text','SHOW');  % DEFAULT
    E.initial_fieldnames = [];
    E.initial_fieldnames = fieldnames(E);
    F.result = E;
    
    % CLUSTER SIZE FILTER %
    E = struct();
    E.item_num = 6;
    E.item_group = 2;
    E.ui_title = uilabel(fig, 'Position',[615 300 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','Cluster Size Filter (px)');
    E.ui_cond = uibuttongroup(fig, 'Position',[530 272 35 40],...
        'SelectionChangedFcn',@(bg,event) recompute_resultmask(fig, E.item_num));
    E.ui_cond_gt = uiradiobutton(E.ui_cond, 'Position',[5 7 30 30], 'Text','>');
    E.ui_cond_lt = uiradiobutton(E.ui_cond, 'Position',[5 -11 30 30], 'Text','<');
    E.ui_cond.SelectedObject = E.ui_cond_gt;  % DEFAULT
    E.mask_cond = E.ui_cond.SelectedObject;
    E.ui_toggle = uibutton(fig, 'state', 'Position',[575 275 35 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_on(fig, E.item_num, btn, false),...
        'Value',false, 'Text','OFF');  % DEFAULT
    E.ui_slider_dec = uibutton(fig, 'push', 'Position',[620 280 20 25], 'Text','-',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, -1, true));
    E.ui_slider = uislider(fig, 'Position',[650 290 200 3],...
        'ValueChangedFcn',@(sld,event) recompute_resultmask(fig, E.item_num, false, false, true),...
        'Limits',[0 100], 'Value',0);  % DEFAULT
    E.ui_slider_inc = uibutton(fig, 'push', 'Position',[860 280 20 25], 'Text','+',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, 1, true));
    E.mask_slider_value = E.ui_slider.Value;
    E.initial_fieldnames = [];
    E.initial_fieldnames = fieldnames(E);
    F.cluster_size = E;
    
    % CLUSTER MERGE DISTANCE %
    E = struct();
    E.item_num = 8;
    E.item_group = 2;
    E.ui_title = uilabel(fig, 'Position',[615 220 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','Cluster Merge Distance (m)');
    E.ui_toggle = uibutton(fig, 'state', 'Position',[575 195 35 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_on(fig, E.item_num, btn, false),...
        'Value',false, 'Text','OFF');  % DEFAULT
    E.ui_slider_dec = uibutton(fig, 'push', 'Position',[620 200 20 25], 'Text','-',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, -1, false));
    E.ui_slider = uislider(fig, 'Position',[650 210 200 3],...
        'ValueChangedFcn',@(sld,event) recompute_resultmask(fig, E.item_num, false, false, true),...
        'Limits',[0 100], 'Value',0);  % DEFAULT
    E.ui_slider_inc = uibutton(fig, 'push', 'Position',[860 200 20 25], 'Text','+',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, 1, false));
    E.mask_slider_value = E.ui_slider.Value;
    E.initial_fieldnames = [];
    E.initial_fieldnames = fieldnames(E);
    F.cluster_dist = E;
    
    % CLUSTER CONCAVITY %
    E = struct();
    E.item_num = 9;
    E.item_group = 2;
    E.ui_title = uilabel(fig, 'Position',[615 140 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','Cluster Concavity');
    E.ui_toggle = uibutton(fig, 'state', 'Position',[575 115 35 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_on(fig, E.item_num, btn, false),...
        'Value',false, 'Text','OFF');  % DEFAULT
    E.ui_slider_dec = uibutton(fig, 'push', 'Position',[620 120 20 25], 'Text','-',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, -1, false));
    E.ui_slider = uislider(fig, 'Position',[650 130 200 3],...
        'ValueChangedFcn',@(sld,event) recompute_resultmask(fig, E.item_num, false, false, true),...
        'Limits',[0.0 1.0], 'Value',0.0);  % DEFAULT
    E.ui_slider_inc = uibutton(fig, 'push', 'Position',[860 120 20 25], 'Text','+',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, 1, false));
    E.mask_slider_value = E.ui_slider.Value;
    E.initial_fieldnames = [];
    E.initial_fieldnames = fieldnames(E);
    F.cluster_concav = E;
    
    % DILATION SIZE %
    E = struct();
    E.item_num = 7;
    E.item_group = 2;
    max_kernel_size = 50;  % SETTING
    E.ui_title = uilabel(fig, 'Position',[615 60 273 20],...
        'VerticalAlignment','Center', 'HorizontalAlignment','Center',...
        'Text','Dilation Kernel Size (sq. side length)');
    E.ui_toggle = uibutton(fig, 'state', 'Position',[575 35 35 35],...
        'ValueChangedFcn',@(btn,event) toggle_item_on(fig, E.item_num, btn, true),...
        'Value',false, 'Text','OFF');  % DEFAULT
    E.ui_slider_dec = uibutton(fig, 'push', 'Position',[620 40 20 25], 'Text','-',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, -1, true, true));
    E.ui_slider = uislider(fig, 'Position',[650 50 200 3],...
        'ValueChangedFcn',@(sld,event) slider_drag_integer(fig, E.item_num),...
        'MajorTicksMode','manual', 'MajorTicks',[0:5:max_kernel_size],...
        'MinorTicksMode','manual', 'MinorTicks',[1:1:max_kernel_size],...
        'Limits',[1 max_kernel_size], 'Value',5);  % DEFAULT
    E.ui_slider_inc = uibutton(fig, 'push', 'Position',[860 40 20 25], 'Text','+',...
        'ButtonPushedFcn',@(btn,event) slider_step(fig, E.item_num, 1, true, true));
    E.mask_slider_value = E.ui_slider.Value;
    if E.ui_toggle.Value == 1
        E.cluster_dilate_size = E.ui_slider.Value;
    else
        E.cluster_dilate_size = 1;
    end
    E.initial_fieldnames = [];
    E.initial_fieldnames = fieldnames(E);
    F.kernel_dilate = E;
    
    guidata(fig, F);
    fig_panel = fig;
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
