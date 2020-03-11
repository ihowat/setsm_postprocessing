function buildSubTiles(tileName,outDir,tileDefFile,databaseFile,waterTileDir,refDemFile)
% buildSubTiles build mosaics from strips in subtiles of 100x100km tiles
%
% buildSubTiles(tileName,outDir,tileDefFile,databaseFile,waterTileDir,refDemFile)
%
% Input arguments:
%tileName=string tile x,y name (e.g. '47_13')
%tileDefFile = matfile list of tile names and ranges (e.g. '/Users/ihowat/unity-home/earthdem/PGC_Imagery_Mosaic_Tiles_Arctic.mat');
%databaseFile = strip database structure produced by by the database function (e.g. '/Users/ihowat/gdrive/projects/earthdem/earthdem_database_unf.mat');
%outDir = directory to ouput subtile files e.g.(['/Users/ihowat/project/howat.4/earthdem/earthdem_mosaic_testing_1km/',tileName]);
%waterTileDir= directory of water/land mask rasters (e.g. '~/data/pgc_projects/ak_water_rasters_v2');
%refDemFile= geotiff reference DEM for quality control ('~/tandemx_alaska_mosaic_3413_tap90m.tif');

% parameters
res=10; % ouput mosaic resolution in meters
subtileSize=1000; % subtile dimensions in meters
buffer=100; % size of tile/subtile boundary buffer in meters
redoFlag = 1; % flag for treating existing subtile .mat files: 1=skip, 2=redo coregistration, 3=redo adjustment
maxNumberOfStrips=100; % maximum number of strips to load in subtile

% if output directory doesnt already exist, make it
if ~exist(outDir,'dir')
    mkdir(outDir)
end

fprintf('Loading tile definition file\n')
tileDefs=load(tileDefFile);
fprintf('Loading strip database and getting tile overlaps\n')
meta=load(databaseFile);

% fix old meta format
if isfield(meta,'f')
    meta.fileName = meta.f;
    meta=rmfield(meta,'f');
end

% get strip areas and alignment stats for quality selection
meta.A = cellfun(@(x,y) polyarea(x,y), meta.x,meta.y);
if isfield(meta,'avg_rmse')
    meta.scene_alignment_meanrmse= meta.avg_rmse;
else
    meta.scene_alignment_meanrmse= cellfun( @(x) mean(x.rmse),...
        meta.scene_alignment);
end

% find index of tile in tile def database
tileInd = find(strcmp(tileDefs.I,tileName));

% get tile boundaries with buffer
x0=tileDefs.x0(tileInd)-buffer;
y0=tileDefs.y0(tileInd)-buffer;
x1=tileDefs.x1(tileInd)+buffer;
y1=tileDefs.y1(tileInd)+buffer;

fprintf('Building water mask tile from %s\n',waterTileDir)
landTile = getTileWaterMask(waterTileDir,tileName,x0,x1,y0,y1,res);

% make array of subtile boundary coordinates
subx0=x0:subtileSize:tileDefs.x1(tileInd)-subtileSize-buffer;
suby0=y0:subtileSize:tileDefs.y1(tileInd)-subtileSize-buffer;
subx1=tileDefs.x0(tileInd)+subtileSize+buffer:subtileSize:x1;
suby1=tileDefs.y0(tileInd)+subtileSize+buffer:subtileSize:y1;

[subx0,suby0] = meshgrid(subx0,suby0);
[subx1,suby1] = meshgrid(subx1,suby1);

subN=numel(subx0);

for n=1:subN
    
    clear x y z offsets za fileNames fa dX dY dZ land
    
    fprintf('subtile %d of %d\n',n,subN)
    
    x=subx0(n):res:subx1(n);
    y=suby1(n):-res:suby0(n);
    y=y(:);
    
    % subset tile land/water mask
    land = interp2(landTile.x,landTile.y(:),single(landTile.z),x,y(:),...
        '*nearest');
    land(isnan(land)) = 0;
    land = logical(land);

    if ~any(land(:))
        fprintf('No land in subtile %d, skipping\n',n)
        continue
    end
    
    outName = [outDir,'/',tileName,'_',num2str(n),'_10m.mat'];
    
    if exist(outName,'file')
        
        % if file exists, make sure it has a za_med array that isnt all
        % NaNs
        mvars  = who('-file',outName);
        
        if any(ismember(mvars,'za_med'))
            
            % yes, has a za_med, load it to check has data
            load(outName,'za_med')
            
            if any(~isnan(za_med(:)))
                
                clear za_med
                
                % yes has data, so treat depending on redoFlag setting
                switch redoFlag
                    case 1 % skip
                        fprintf('%s exists, skipping\n',outName)
                        continue
                    case 2 % redo starting at coregistratiom
                        fprintf('%s exists, restarting at coregistration\n',outName)
                        load(outName,'x','y','z','t','fileNames');
                    case 3 % redo starting at adjustment
                        fprintf('%s exists,  restarting at adjustment\n',outName)
                        load(outName,'x','y','z','t','fileNames','offsets');
                    case 4 % redo at output
                        fprintf('%s exists,  restarting at median ouput\n',outName)
                        load(outName,'x','y','za','t','fileNames','dZ','dX','dY','fa');
                        % dummy vars to skip check points
                        z = [];
                        offsets=[];
                    case 5 % write 2m
                        fprintf('%s exists, re-writing 2m output\n',outName)
                        load(outName)
                        if ~exist('fileNames0','var')
                            fileNames0 = fileNames;
                        end
                end
            else
                load(outName)
                clear za za_med dZ dX dY
            end
            
        else
            load(outName)
            clear za za_med dZ dX dY
            
        end
    end
    
    % find overlapping strips
    if ~exist('fileNames','var')
        
        % make sure this subtile has land surface in it
        subtilePoly = polyshape([subx0(n);subx0(n);subx1(n);subx1(n)],...
            [suby0(n);suby1(n);suby1(n);suby0(n)]);
        
        % index of strips in this subtile
        ind=stripSearch(meta.x,meta.y,subtilePoly);
        
        % check for maximum #'s of overlaps
        if length(ind) > maxNumberOfStrips
            
            %first remove single subscenes
            ind(meta.scene_alignment_meanrmse(ind) == 0 |...
                isnan(meta.scene_alignment_meanrmse(ind))) = [];
            
        end
        
        % check if still more strips than max
        if length(ind) >  maxNumberOfStrips
            
            % sort by weigthed area and coverage
            scene_alignment_meanrmse_scaled= -(meta.scene_alignment_meanrmse(ind)-...
                mean(meta.scene_alignment_meanrmse(ind)))./...
                std(meta.scene_alignment_meanrmse(ind));
            
            % normalize areas for group
            A_scaled = (meta.A(ind)-mean(meta.A(ind)))./std(meta.A(ind));
            
            % normalize scene alignment stats for group
            scene_alignment_meanrmse_scaled = scene_alignment_meanrmse_scaled +...
                abs(min(scene_alignment_meanrmse_scaled));
            
            % take the sum of normalized values
            A_scaled = A_scaled + abs(min(A_scaled));
            
            % sort by scaled sum
            [~,ind_rank] = sort(A_scaled.*scene_alignment_meanrmse_scaled,...
                'descend');
            
            % take top N strips in sort
            ind = ind(ind_rank(1:maxNumberOfStrips));
        end
        
        if isempty(ind)
            fprintf('no overlapping strips, skipping\n')
            continue
        end
        
        % convert meta names into dem names
        fileNames  = strrep(meta.fileName(ind),'meta.txt','dem_10m.tif');
        
        % convert filenames for mounted drives on laptop
        if ismac
            fileNames  = strrep(fileNames,'/fs/project/howat.4',...
                '/Users/ihowat/project');
            fileNames  = strrep(fileNames,'/data',...
                '/Users/ihowat/data');
        else
            fileNames  = strrep(fileNames,'/data',...
                '/home/howat.4/data');
            fileNames = strrep(fileNames,'/mnt/pgc/home/howat.4/data',...
                '/mnt/pgc/data');
        end
        
        fileNames = strrep(fileNames,'meta.txt','dem_10m.tif');
        
        if  ~exist(fileNames{1},'file')
            
              fileNames = strrep(fileNames,'2m_dem_10m.tif','10m_dem.tif');
                 
                  if  ~exist(fileNames{1},'file')
                      error('1 or more files dont exist: %s\n', fileNames{1});
                  end
        end
        
        % make date vector
        [~,name] =  cellfun(@fileparts,fileNames,'uniformoutput',0);
        t=cellfun(@(x) datenum(x(6:13),'yyyymmdd'),name)';
        
    end
    
    % extract strip subsets into stack
    if ~exist('z','var')
        
        fprintf('extracting %d strip subsets ',length(fileNames))
        [x,y,z,missingFlag] =extractSubGrid(fileNames,subx0(n),subx1(n),...
            suby0(n),suby1(n),res);
        
        % remove layers missing data
        z = z(:,:,~missingFlag);
        fileNames=fileNames(~missingFlag);

        % save full filename list with repeat segments for 2m
        fileNames0 = fileNames;
        
        % merge segmentsfrom same strips
        % dont get offsets between segs in same strip:make a vector of z's
        % that ar belonging to the same strip
        [~,stripid] =  cellfun(@fileparts,fileNames,'uniformoutput',0);
        stripid =  cellfun(@(x) x(1:47),stripid,'uniformoutput',0);
        unique_stripid = unique(stripid);
        
        r=[];
        it=1;
        for it=1:length(unique_stripid)
            
            segs = find(strcmp(stripid,unique_stripid(it)));
            
            segs=segs(:); 
            
            if length(segs) > 1
                
                z(:,:,segs(1)) = nanmedian(z(:,:,segs),3);

                r = [r(:);segs(2:end)];
            end
        end
        
        z(:,:,r) = [];
        fileNames(r) = [];
        t(r) = [];
        
        clear it segs r
        
        fprintf('saving x,y,z,t,fileNames to %s\n',outName)
        save(outName,'x','y','z','t','fileNames','fileNames0','-v7.3');
        
    end
    
    % calculate DEM pairwise offsets
    if ~exist('offsets','var')
        
        if ~exist('land','var')
            % subset tile land/water mask
            land = interp2(landTile.x,landTile.y(:),single(landTile.z),x,y(:),...
                '*nearest');
            land(isnan(land)) = 0;
            land = logical(land);
        end
        
        if ~any(land(:))
            fprintf('No land in subtile %d, skipping\n',n)
            continue
        end
        
        fprintf('saving land to %s\n',outName)
        save(outName,'land','-append');
        
        % dont get offsets between segs in same strip:make a vector of z's
        % that ar belonging to the same strip
        [~,stripid] =  cellfun(@fileparts,fileNames,'uniformoutput',0);
        stripid =  cellfun(@(x) x(1:47),stripid,'uniformoutput',0);
        [~,~,strip_ind] = unique(stripid);
        
        fprintf('performing pairwise coregistration, ')
        offsets=coregisterStack(x,y,z,land,strip_ind);
        
        fprintf('saving offsets to %s\n',outName)
        save(outName,'offsets','-append');
        
    end
    
    % adjustment
    if ~exist('dZ','var')
        
        % set offsets with horizontal shift failure to zero
        offsets.dx(offsets.dxe == 0) = NaN;
        offsets.dy(offsets.dye == 0) = NaN;
        
        offsets.dxe(isnan(offsets.dx)) = NaN;
        offsets.dye(isnan(offsets.dy)) = NaN;
        
        
        if ~any(~isnan(offsets.dz))
            
            % if all offsets failed, set adjustment thresholds to inf and
            % shifts to zero
            offsetErrMax = inf;
            min_sigma_dz_coregMax= inf;
            min_abs_mean_dz_coregMax= inf;
            min_abs_median_dz_coregMax = inf;
            
            dZ = zeros(size(z,3),1);
            dX = zeros(size(z,3),1);
            dY = zeros(size(z,3),1);
            
        else
            
            %number of grid points with coverage
            c0 = any(z,3);
            c0 = sum(c0(:));
            
            % iteratively perform adjustment, relaxing thresholds by 10%, up to a
            % maximum of 500%  , until coverage is maximum
            it=1;
            while it < 5
                
                % pairwise coregistration statistics filter threshold defaults
                offsetErrMax = 0.1*it;
                min_sigma_dz_coregMax=4*it;
                min_abs_mean_dz_coregMax=0.1*it;
                min_abs_median_dz_coregMax = 1*it;
                
                % perform adjustment for this iteration
                [dZ,dX,dY] = adjustOffsets(offsets,'offsetErrMax',offsetErrMax,...
                    'min_sigma_dz_coregMax',min_sigma_dz_coregMax,...
                    'min_abs_mean_dz_coregMax',min_abs_mean_dz_coregMax,...
                    'min_abs_median_dz_coregMax',min_abs_median_dz_coregMax);
                
                %check for spatial coverage of z layers with solutions
                c1 = any(z(:,:,~isnan(dZ)),3);
                c1 = sum(c1(:));
                
                % if coverage is 100%, break
                if c1 == c0
                    break
                end
                
                % if coverage is < 100%, increase thresholds by 10% and
                % repeat
                it = it+0.1;
            end
            
            % If adjustment fails, take single DEM with best comparison to
            % reference DEM
            if ~any(~isnan(dZ))
                
                % load subset of reference dem covering subtile
                zr=readGeotiff(refDemFile,...
                    'map_subset',[x(1)-90 x(end)+90 y(end)-90 y(1)+90]);
                
                % interpolate reference dem to subtile
                try
                    zri = interp2(zr.x,zr.y(:),zr.z,x,y(:),'*linear');
                catch
                    fprintf('Reference dem failed to interpolate to subtile, skipping\n');
                    continue
                end
                
                % create vertical std dev vector
                dz_std = nan(size(z,3),1);
                
                % loop through dems in subtile stack and calculate vertical
                % difference between dem and reference dem, saving the
                % standard deviation over land.
                for i=1:size(z,3)
                    dz =  zri - z(:,:,i);
                    dz_std(i) = nanstd(dz(land));
                end
                
                % sort the dems by lowest std dev from reference and select
                % the least number of dems needed for 100% cover
                [~,nsort] = sort(dz_std,'ascend');
                
                % successively calculate coverage provided by each ith dem
                % in order of increasing std dev from reference, breaking
                % when coverage is 100%
                for i=1:length(nsort)
                    c = any(z(:,:,nsort(1:i)),3);
                    if sum(c(:)) == c0
                        break
                    end
                end
                
                % only retain the top ith dems
                nsort = sort(nsort(1:i));
                
                % set adjustment thresholds to inf to shut them off or
                % indicate no adjustment
                offsetErrMax = inf;
                min_sigma_dz_coregMax= inf;
                min_abs_mean_dz_coregMax= inf;
                min_abs_median_dz_coregMax = inf;
                
                % only perform adjustment if more than 1 dem
                if length(nsort) > 1
                    
                    % find coregistration offsets for pairs of nsort dems
                    in = ismember(offsets.i,nsort) & ismember(offsets.j,nsort);
                    
                    % if these DEMs dont overlap, or are of the same strip,
                    % they wont have offsets: skip
                    if any(in)
                        
                        % sample the offset structure for the pairs of nsort
                        offsets_sub = structfun(@(x) x(in), offsets,'uniformoutput',0);
                        
                        % perform adjustment for pairs of nsort
                        [dZ(nsort),dX(nsort),dY(nsort)] = adjustOffsets(offsets_sub,...
                            'offsetErrMax',offsetErrMax,...
                            'min_sigma_dz_coregMax',min_sigma_dz_coregMax,...
                            'min_abs_mean_dz_coregMax',min_abs_mean_dz_coregMax,...
                            'min_abs_median_dz_coregMax',min_abs_median_dz_coregMax);
                    end
                end
                
                % if adjustment still fails, or just one dem, just set to
                % zero for nsort "best" dems.
                if ~any(~isnan(dZ))
                    dZ(nsort) = 0;
                    dX(nsort) = 0;
                    dY(nsort) = 0;
                end
                
            end
            
        end
        
        fprintf('saving adjustment variables to %s\n',outName)
        save(outName,'dZ','dX','dY','offsetErrMax','min_sigma_dz_coregMax',...
            'min_abs_mean_dz_coregMax','min_abs_median_dz_coregMax','-append');
        
    end
    
    % apply adustment
    if ~exist('za','var')
        
        % layers with missing adjustments
        n_missing = isnan(dZ);
        
        % make adjusted z array
        za=nan(size(z),'single');
        
        % index of layers with adjustments
        iterVec = 1:size(z,3);
        iterVec(n_missing) = [];
        
        fprintf('applying adjustment\n')
        clear k
        for k=iterVec
            za(:,:,k) = interp2(x + dX(k),y + dY(k), z(:,:,k) + dZ(k),...
                x,y,'*linear');
        end
        
        fprintf('saving za to %s\n',outName)
        save(outName,'za','-append');
        
    end
    
    if ~exist('fa','var')
        % apply a pixel-by-pixel filter to remove outliers
        fa = pairwiseDifferenceFilter(za,'mask',land,'minmad',2);
        fprintf('saving fa to %s\n',outName)
        save(outName,'fa','-append');
    end
    
    % apply filter and get median
    za(~fa) = NaN;
    
    za_med = single(nanmedian(za,3));
    N = uint8(sum(~isnan(za),3));
    
    fprintf('saving za_med and N to %s\n',outName)
    save(outName,'za_med','N','-append');
    
    fprintf('making 2m version\n')
    outName2m = strrep(outName,'_10m.mat','_2m.mat');
    
    % if strip segments were combined, need to expand offset vectors and fa
    % array to match orginal file list
    if length(fileNames0) ~= length(fileNames) 

        [~,stripid] =  cellfun(@fileparts,fileNames0,'uniformoutput',0);
        stripid =  cellfun(@(x) x(1:47),stripid,'uniformoutput',0);
        [~,~,strip_ind] = unique(stripid);
         
         dZ = dZ(strip_ind);
         dX = dX(strip_ind);
         dY = dY(strip_ind);
         
         fa = fa(:,:,strip_ind);
    
    end
    
    make2m(fileNames0,x,y,dZ,dX,dY,land,fa,outName2m);
    
end

function make2m(fileNames,x,y,dZ,dX,dY,land,fa,outName)

% make date vector
[~,name] =  cellfun(@fileparts,fileNames,'uniformoutput',0);
t=cellfun(@(x) datenum(x(6:13),'yyyymmdd'),name)';

% layers with missing adjustments
n_missing = isnan(dZ);

dZ(n_missing) = [];
dX(n_missing) = [];
dY(n_missing) = [];
fileNames(n_missing) = [];
fa(:,:,n_missing) = [];
t(n_missing) = [];

fileNames = strrep(fileNames,'_10m.tif','.tif');

[x,y,z,~,mt] =extractSubGrid(fileNames,min(x),max(x),...
    min(y),max(y),2);

% merge segmentsfrom same strips
% dont get offsets between segs in same strip:make a vector of z's
% that ar belonging to the same strip
[~,stripid] =  cellfun(@fileparts,fileNames,'uniformoutput',0);
stripid =  cellfun(@(x) x(1:47),stripid,'uniformoutput',0);
unique_stripid = unique(stripid);

r=[];
it=1;
for it=1:length(unique_stripid)
    
    segs = find(strcmp(stripid,unique_stripid(it)));
    
    if length(segs) > 1
        
        z(:,:,segs(1)) = nanmedian(z(:,:,segs),3);
        
        mt(:,:,segs(1)) = any(mt(:,:,segs),3);
        
        r = [r,segs(2:end)];
    end
end

z(:,:,r) = [];
mt(:,:,r) = [];
fileNames(r) = [];
t(r) = [];
dZ(r) = [];
dX(r) = [];
dY(r) = [];

clear it segs r

% make adjusted z and mt arrays
za=nan(size(z),'single');
mta = false(size(mt));
for k=1:size(z,3)
    zak = interp2(x + dX(k),y + dY(k), z(:,:,k) + dZ(k),...
        x,y,'*linear');
    
    mtak = interp2(x + dX(k),y + dY(k), single(mt(:,:,k)),...
        x,y,'*nearest');
    
    mtak(isnan(mtak)) = 0;
    mtak = logical(mtak);
    
    if any(any(~fa(:,:,k)))
        fak = imresize(fa(:,:,k),size(za(:,:,k)),'nearest');
        zak(~fak)=NaN;
        mtak(~fak) = false;
    end
    
    za(:,:,k) = zak;
    mta(:,:,k) = mtak;
end

za_med = nanmedian(za,3);
%za_std =  nanstd(za,[],3);
za_mad = mad(za,1,3);
N = uint8(sum(~isnan(za),3));
Nmt = uint8(sum(mta,3));

% resize land mask
land = imresize(land,size(za_med),'nearest');

t=t-datenum('1/1/2000 00:00:00');
t=reshape(t,1,1,[]);
t=repmat(t,size(za_med));
t(isnan(za))=NaN;
tmax = max(t,[],3);
tmin = min(t,[],3);
%tmean = nanmean(t,3);

tmax = uint16(tmax);
tmin = uint16(tmin);
%tmean = uint16(tmean);


% Incomplete attempt at code for retrieving dates of median values
% [za_sort,n]  = sort(za,3);
% isodd=logical(mod(single(N),2));
% n1=zeros(size(N),'uint8');
% n2=zeros(size(N),'uint8');
%
% n1(isodd & N > 0) = uint8(ceil(single(N(isodd & N > 0))./2));
% n2(isodd & N > 0) = n1(isodd & N > 0);
%
% n1(~isodd & N > 0) = uint8(single(N(~isodd & N > 0))./2);
% n2(~isodd & N > 0) = n1(~isodd & N > 0)+1;
%
% [col,row] = meshgrid((1:size(za,2))',1:size(za,1));
%
% n1(N == 0) = [];
% n2(N == 0) = [];
% row(N == 0) = [];
% col(N == 0) = [];
%
% ind1 = sub2ind(size(za),row(:),col(:),n1(:));
% ind2 = sub2ind(size(za),row(:),col(:),n2(:));
%
% za1 = za_sort(ind1);
% za2 = za_sort(ind2);
%
%
% za_med = nan(size(za,1),size(za,2),'single');
% ind = sub2ind(size(za_med),row(:),col(:));
% za_med(ind) = (za1+za2)./2;
%
% tmed= zeros(size(za,1),size(za,2),'uint16');
% t_med(ind) = (ta1+ta2)./2;
%
% t =
%
% za1 = za_sort(ind1);
% za2 = za_sort(ind2);

fprintf('saving x, y za_med land za_mad N tmax tmin to %s\n',outName)
save(outName,'x','y','za_med','land','za_mad','N','Nmt','tmax','tmin','-v7.3');

