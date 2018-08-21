function qcTile(tileFile,meta)
% qcTile: quality control REMA/ArcticDEM tile mosaic
% qcTile(tileFile,meta) where tileFile contains the file path/name of the
% tile and meta is the meta data structure. The tileFile must end with
% *dem.mat. If a file named *dem_shade.tif exists, that will be used for
% the qc display and, if not, a hillshade will be created. Hillshades must
% exist for the individual strips, however.

global prompt
global qc
global IA
global poly_handles
global DELETE_ENABLED
global QC_LOCK_FILE
global USERNAME
DELETE_ENABLED = false;
QC_LOCK_FILE = '';
USERNAME = getenv('username');
LOCK_FILE_CLEANUP = onCleanup(@() remove_lock());

tileschema='V:/pgc/data/scratch/claire/repos/setsm_postprocessing_pgc/PGC_Imagery_Mosaic_Tiles_Arctic.mat'; %PGC/NGA Tile definition file, required
changePath= 'V:/pgc'; %if set, will change the path to the REMA directory from what's in the database file. set to [] if none.

dbasedir_local = [getenv('USERPROFILE'),'\setsm_postprocessing_dbase'];
[tileschema] = copy_dbase_local(dbasedir_local, tileschema);

fprintf('loading tile schema\n');
tiles=load(tileschema);

% get tile with matching name
[~,tn,~] = fileparts(tileFile);
tn=strrep(tn,'_8m_dem','');
i=strcmp(tn,tiles.I);
tile = structfun(@(x) ( x(i) ), tiles, 'UniformOutput', false);

% open the tileFile into a matfile structure
m=matfile(tileFile);

% check if shade tif exists and read if it does
% fprintf('loading/building tile hillshade\n');
shadeFile = strrep(tileFile,'dem.mat','dem_shade.tif');
if strcmp(tileFile,shadeFile); shadeFile =[]; end
if exist(shadeFile,'file')
    fprintf('loading tile hillshade\n');
    h = readGeotiff(shadeFile);
    x = h.x;
    y = h.y;
    h = h.z;
    
else
    % no shade so create one
    fprintf('building tile hillshade\n');
    h=hillshade(m.z,m.x,m.y);
    x = m.x;
    y = m.y;
end

% plot the hillshade in the current figure window
f1 = gcf;
imagesc(x,y,uint8(h));
axis xy equal;
colormap gray
set(gca,'Position',[0 0 1 1]);
set(gca,'clim', [0 255])
hold on

if isfield(tile,'coastline')
    i=1;
    for i=1:length(tile.coastline{1})
        if isempty(tile.coastline{1}{i}); continue; end
        plot(tile.coastline{1}{i}(1,:),...
            tile.coastline{1}{i}(2,:),'b','linewidth',1)
    end
end

% build coverage grid
res = 40;
buff = 0;

coverage_tile_x = tile.x0-buff*res: res:tile.x1+buff*res;
coverage_tile_y = tile.y1+buff*res:-res:tile.y0-buff*res;
coverage_tile_y = coverage_tile_y(:);


deleteFlag = false;
while true
    
    set(0, 'currentfigure', f1);
    
    qc_flag_filter = [];
    while true
        prompt = 'Draw ROI? (y/n/[1-5] qc flag filter)\n';
        s=input(prompt,'s');

        if strcmpi(s,'y') || strcmpi(s,'n') || strcmpi(s,'.') || strcmpi(s,'0')
            break;
        elseif ~isempty(s) && all(ismember(s, '12345'))
            qc_flag_filter = str2num(s);
            break;
        else
            fprintf('%s not recogized, try again\n',s);
        end
    end
    
    if strcmpi(s,'n') || strcmpi(s,'0')
        break;
        
    elseif strcmpi(s,'c')
         clim(1) =input('enter min c value\n');
         clim(2) = input('enter max c value\n');
        
         set(gca,'clim',clim);
         
    else
        
        deleteFlag = true;
        
        if ~exist('f','var')
            %% match strip files to database
            f=m.f; % extract image filename vector from from matfile
%             f(isnan(m.rmse))=[]; % remove unused data
            [~,IA]= intersect(meta.f,f); % find index of tile files in meta
            meta = structfun(@(x) ( x(IA,:) ), meta, 'UniformOutput', false); %crop meta
            meta = addQC2Meta(meta); % add existing qc data
        end
        
        [BW_roi,xv,yv] =  roipoly; % draw ROI
        if isempty(BW_roi)
            continue;
        end
        BW_roi=interp2(x,y,BW_roi,coverage_tile_x,coverage_tile_y,'*nearest');
        plot(xv,yv,'b','linewidth',2); % plot ROI
        
        % find strips overlapping ROI
        n = find(meta.xmax > min(xv) &  meta.xmin < max(xv) & ...
            meta.ymax > min(yv) &  meta.ymin < max(yv));
        
        in=zeros(size(n));
        for i=1:length(in)
            in(i) = any(inpolygon(meta.x{n(i)},meta.y{n(i)},xv,yv)) | ...
                any(inpolygon(xv,yv,meta.x{n(i)},meta.y{n(i)}));
        end
        n=n(logical(in));
        
        % select the whichever registration has the better sigma_bias (all or 1 yr)
        if isfield(meta,'sigma_all') &&  isfield(meta,'sigma_1yr') &&  ~isfield(meta,'sigma')
            meta.sigma = nanmin([meta.sigma_all(:)';meta.sigma_1yr(:)'])';
            meta.sigma(meta.sigma > 1) = NaN;
        else
            meta.sigma  = nan(size(meta.f));
        end

        meta.avg_rmse(meta.avg_rmse == 0) = NaN;
        
        
        % order strips that overlap roi by number of data pixels in roi
        overlap_size = zeros(size(n));
        
        % can be slow, so we'll use a counter
        count = 0;
        
        % file loop
        for v=1:length(n)
            i = n(v);

            % counter
            if v>1; for p=1:count; fprintf('\b'); end; end %delete line before
            count = fprintf('strip %d of %d',v,length(n));

            % load strip data coverage
            coverageFile = strrep(meta.f{i}, 'meta.txt', 'dem_coverage.tif');
            if ~isempty(changePath)
                coverageFile=strrep(coverageFile,'/mnt/pgc',changePath);
                coverageFile=strrep(coverageFile,'/','\');
            end
            if exist(coverageFile, 'file')
                
                BW_coverage = get_tile_size_coverage(coverageFile, coverage_tile_x, coverage_tile_y);
                % if mask data exists, apply it
                if meta.qc(i) == 3
                    for j=1:length(meta.maskPolyx{i})
                        BW_coverage(roipoly(x,y,BW_coverage,...
                            meta.maskPolyx{i}{j},meta.maskPolyy{i}{j}))=0;
                    end
                end
                
                overlap_size(v) = nnz(BW_roi & BW_coverage);
                
            else
                overlap_size(v) = Inf;
            end
            
            % get rid of mask
            clear BW_coverage
        end
        
        % get rid of mask
        clear BW_roi
        
        [strip_data_pixels_overlapping_roi,overlap_sort_I] = sort(overlap_size, 'descend');
        n = n(overlap_sort_I);
        strip_data_pixels_overlapping_roi
        
        
        f2=figure;
        F = [];
        i=0;
        while true
            i=i+1;
            if i > length(n)
                break;
            end
            
            set(0, 'currentfigure', f2);
            clf;
            
            poly_handles = [];
            DELETE_ENABLED = false;
     
            % load and plot strip hillshade
            fileName=strrep(meta.f{n(i)},'meta.txt','dem_browse.tif');
            if ~isempty(changePath)
                fileName=strrep(fileName,'/mnt/pgc',changePath);
                fileName=strrep(fileName,'/','\');
            end
            orthoFile=strrep(fileName,'dem_browse.tif','ortho_browse.tif');
            fprintf('loading %d of %d: %s\n',i,length(n),fileName);
            
            % load qc data
            qcFile = [fileparts(fileName),'\qc.mat'];
            qcLockFile = [qcFile,'.lock'];
            if exist(qcLockFile, 'file') ~= 2
                fprintf(2, 'Please have an admin generate an empty qc.mat.lock file for this region with proper ACL setting: %s\n', qcLockFile);
            else
                user_id = read_lock(qcLockFile);
                if isempty(user_id)
                    remove_lock();
                    QC_LOCK_FILE = qcLockFile;
                    write_lock();
                elseif strcmp(user_id, USERNAME) && strcmp(qcLockFile, QC_LOCK_FILE)
                    ;
                else
                    fprintf(2, 'lock file exists for user "%s", skipping: %s\n', user_id, qcLockFile);
                    fprintf(2, 'if this region is locked in error, simply OPEN the lock file in a text editor, CLEAR its contents, and SAVE to unlock\n');
                    continue;
                end
            end
            qc=load(qcFile);
            fileNames=strrep(qc.fileNames,'/','\');
            [~,IA]=intersect(fileNames, fileName);
            
            if ~isempty(qc_flag_filter) && qc.flag(IA) ~= qc_flag_filter
                continue;
            end
            
            I=readGeotiff(fileName);
            
            
            if isempty(F)
                F = struct();
            else
                F = guidata(f2);
                gui_fields = fieldnames(F);
                for k = 1:length(gui_fields)
                    F = rmfield(F, gui_fields{k});
                end
            end
            F.plot_group_hill = hggroup;
            F.plot_group_ortho = hggroup;
            F.fileName_ortho = orthoFile;
            F.ortho_browse_ready = false;
            guidata(f2, F);

            % Make square image to show in the figure window.
            Z_fig = get_square_browse_images(f2, I);
            F = guidata(f2);
            
            
            imagesc(F.X_fig,F.Y_fig,Z_fig,'alphadata',single(Z_fig ~= 0), 'Parent',F.plot_group_hill)
            set(gca,'color','r')
            set(gca,'Position',[0 0 1 1]);
            axis xy  equal tight;
            colormap gray;
%             set(gcf,'units','normalized');
%             set(gcf,'position',[0.01,0.01,.35,.9])
            hold on
            % plot ROI
            plot(xv,yv,'b','linewidth',2)
            set(gca,'xlim',[min(F.X_fig)-500 max(F.X_fig)+500],'ylim',[min(F.Y_fig)-500 max(F.Y_fig)+500])
            
            
            % Make radio buttons for switching between hillshade, ortho.
            F.ui_image_select = uibuttongroup('Visible','off',...
                              'Units','pixels',...
                              'Position',[0 0 130 100],...
                              'SelectionChangedFcn',@(source,event) change_browse_image(f2));
            uicontrol(F.ui_image_select,'Style','radiobutton', 'Position',[20 50 100 30],...
                              'String','DEM Hillshade', 'HandleVisibility','off');
            uicontrol(F.ui_image_select,'Style','radiobutton', 'Position',[20 20 100 30],...
                              'String','Ortho Image', 'HandleVisibility','off');
            F.ui_image_select.Visible = 'on';
            guidata(f2, F);
            
            
            % existing qc is 3, plot the mask polys
            j=1;
            if qc.flag(IA) == 3
                while true
                    if j > length(qc.x{IA})
                        break;
                    elseif isempty(qc.x{IA}{j})
                        qc.x{IA}(j) = [];
                        qc.y{IA}(j) = [];
                        continue;
                    end
                    hl = plot(qc.x{IA}{j},qc.y{IA}{j},'g','linewidth',2);
                    set(hl, 'ButtonDownFcn', @(~,~) delete_poly(j), 'HitTest','on');
                    poly_handles = [poly_handles, hl];
                    j=j+1;
                end
            end
            
            
            %% specify flag
            j=1;
            while j
                
                fprintf('gcp sigma=%.2f, mean coreg RMSE=%.2f, max coreg RMSE=%.2f, current qc flag: %d\n',...
                    meta.sigma(n(i)),meta.avg_rmse(n(i)), meta.max_rmse(n(i)), qc.flag(IA));
                
                prompt = 'Enter quality flag: ENTER=NO CHANGE, 0=BACK, 9=EXIT ROI, 1=good, 2=partial, 3=manual edit, 4=poor, 5=unuseable\n';
                s=input(prompt,'s');
                
                if isempty(s); break; end
                
                if ~isempty(s) && all(ismember(s, '0123456789'))
                    flag = str2num(s);
                    if flag == 1 || flag == 2 || flag == 3 || flag == 4 || flag == 5 || flag == 0 || flag == 9
                        break;
                    end
                end
                
                fprintf('%s not recogized, try again\n',s);
                
            end
            
            if isempty(s); continue; end
            
            if flag == 0
                i = max(i-2, 0);
                continue;
            end
            
            if flag == 9
                break;
            end
            
            %% if manual selected 
            if flag ==3
                  
                
                j=1;
                if qc.flag(IA) == 3
                    
%                     prompt = 'keep existing mask polygons? (y/n)\n';
%                     s=input(prompt,'s');
%                     
%                     if strcmpi(s,'n')
%                         
%                         qc.x{IA}=cell(1); qc.y{IA}=cell(1);
%                         delete(poly_handles);
%                         poly_handles = [];
%                     else
%                         
%                         j= length(qc.x{IA})+1;
%                     end
                    
                    j= length(qc.x{IA})+1;
                    
                end

                fprintf('entering manual edit mode\n');
                fprintf(2, '-- clicking on a polygon edge will prompt for deletion --\n');
                DELETE_ENABLED = true;
                
                while j
                    
                    while j
                        prompt = 'continue editing this image? (y/n)\n';
                        s=input(prompt,'s');
                        
                        if ~strcmpi(s,'y') && ~strcmpi(s,'n') && ~strcmpi(s,'.') && ~strcmpi(s,'0')
                            fprintf('%s not recogized, try again\n',s);
                        else
                            break
                        end
                    end
                    
                    if strcmpi(s,'n') || strcmpi(s,'0'); break; end
                    
                    mask_region_num = 1;
                    for k=1:length(qc.x{IA})
                        if ~isempty(qc.x{IA}{k})
                            mask_region_num = mask_region_num + 1;
                        end
                    end
                    fprintf('drawing mask region %d\n',mask_region_num);
                    
                     
                    figure(f2);
                    [roi_BW,roi_x,roi_y] = roipoly;
                    if isempty(roi_BW)
                        continue;
                    end
                    qc.x{IA}{j} = roi_x;
                    qc.y{IA}{j} = roi_y;
                    clear roi_BW;

                    hl = plot(qc.x{IA}{j},qc.y{IA}{j},'g','linewidth',2);
                    set(hl, 'ButtonDownFcn', @(~,~) delete_poly(length(qc.x{IA})), 'HitTest','on');
                    poly_handles = [poly_handles, hl];
                    
                    j=j+1;
                 end      
            end
            
            
            qc.flag(IA)=flag;
            
            if flag == 1 || flag == 2
                qc.x{IA}=cell(1); qc.y{IA}=cell(1);
                
            end
            
            if flag == 4 || flag == 5
                qc.x{IA}=cell(1); qc.y{IA}=cell(1);
                qc.x{IA}{1} = [min(I.x);min(I.x);max(I.x);max(I.x);min(I.x)];
                qc.y{IA}{1} = [min(I.y);max(I.y);max(I.y);min(I.y);min(I.y)];
            end
            
            
            save(qcFile,'-struct','qc');
            
        end
        
        close(f2);
              
    end
    
end

% if deleteFlag
% 
%     s=input('Delete edited tile file (y/n)\n','s');
%     
%     
%     if strcmpi(s,'y')
%         
%          [filePath,fileName] = fileparts(tileFile);
%          delete([filePath,'/',fileName(1:5),'*']) 
%     end
% end

clf
end


function delete_poly(index_num)
    global prompt
    global qc
    global IA
    global poly_handles
    global DELETE_ENABLED
    
    if ~DELETE_ENABLED
        return;
    end
    
    for k=1:length(qc.x{IA})
        set(poly_handles(k), 'ButtonDownFcn', []);
    end
    
    hl = plot(qc.x{IA}{index_num},qc.y{IA}{index_num},'w','linewidth',3);
    
    while true
        fprintf('delete highlighted polygon? (y/n)\n');
        s=input('','s');

        if ~strcmpi(s,'y') && ~strcmpi(s,'n') && ~strcmpi(s,'.') && ~strcmpi(s,'0')
            fprintf('%s not recogized, try again\n',s);
        else
            break
        end
    end
    
    if strcmpi(s,'y') || strcmpi(s,'.')
        qc.x{IA}{index_num} = [];
        qc.y{IA}{index_num} = [];
        
        set(poly_handles(index_num), 'visible','off');
        
        fprintf('polygon removed\n');
    end
    
    delete(hl);
    
    for k=1:length(qc.x{IA})
        set(poly_handles(k), 'ButtonDownFcn', @(~,~) delete_poly(k), 'HitTest','on');
    end
    
    fprintf(prompt);

end


function write_lock()
    global QC_LOCK_FILE
    global USERNAME
    
    fid = fopen(QC_LOCK_FILE, 'wt');
    fwrite(fid, USERNAME);
    fclose(fid);

end


function [user_id] = read_lock(qcLockFile)

    lock_txt = fileread(qcLockFile);
    user_id = strrep(lock_txt, char(10), '');
    user_id = strrep(user_id, '\n', '');
    
end


function remove_lock()
    global QC_LOCK_FILE
    
    if exist(QC_LOCK_FILE, 'file') == 2
        fid = fopen(QC_LOCK_FILE, 'w');
        fclose(fid);
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


function [Z_fig] = get_square_browse_images(fig, I)
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
    
    Z_fig = zeros(max(size(Y_fig)), max(size(X_fig)));
    
    Z_c0 = find(I.x(1) == X_fig);
    Z_c1 = find(I.x(end) == X_fig);
    Z_r0 = find(I.y(1) == Y_fig);
    Z_r1 = find(I.y(end) == Y_fig);
    
    Z_fig(Z_r0:Z_r1, Z_c0:Z_c1) = I.z;
    
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


function [dstfile] = make_local_copy(srcfile, dstdir)
    srcfile_stats = dir(srcfile);
    srcfname = srcfile_stats.name;
    dstfile = [dstdir,'\',srcfname];
    
    if exist(dstdir, 'dir') ~= 7
        fprintf('Making local copy of db file at %s ...', dstfile);
        mkdir(dstdir);
        
    elseif exist(dstfile, 'file') == 2
        dstfile_stats = dir(dstfile);
        if dstfile_stats.datenum == srcfile_stats.datenum
            % Local copy doesn't needed to be updated.
            return;
        else
            fprintf('Updating local copy of db file at %s ...', dstfile);
        end
        
    else
        fprintf('Making local copy of db file at %s ...', dstfile);
    end
    
    status = copyfile(srcfile, dstdir);
    if status == 0
        fprintf(' failed!\n');
        dstfile = '';
    elseif status == 1
        fprintf(' success!\n');
    end
end


function [varargout] = copy_dbase_local(dbasedir_local, varargin)
    for i = 1:length(varargin)
        dbasefile = varargin{i};
        dbasefile_local = make_local_copy(dbasefile, dbasedir_local);
        if isempty(dbasefile_local)
            fprintf('Falling back to network load\n');
        else
            dbasefile = dbasefile_local;
        end
        varargout{i} = dbasefile;
    end
end
