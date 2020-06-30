function buildSubTilesSortFix(tileName,outDir,tileDefFile,databaseFile,waterTileDir,refDemFile)
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
%
% % parameters
% res=10; % ouput mosaic resolution in meters
% subtileSize=1000; % subtile dimensions in meters
% buffer=100; % size of tile/subtile boundary buffer in meters
% redoFlag = 1; % flag for treating existing subtile .mat files: 1=skip, 2=redo coregistration, 3=redo adjustment
% maxNumberOfStrips=100; % maximum number of strips to load in subtile
%
% % if output directory doesnt already exist, make it
% if ~exist(outDir,'dir')
%     mkdir(outDir)
% end
%
% fprintf('Loading tile definition file\n')
% tileDefs=load(tileDefFile);
% fprintf('Loading strip database and getting tile overlaps\n')
% meta=load(databaseFile);
%
% % fix old meta format
% if isfield(meta,'f')
%     meta.fileName = meta.f;
%     meta=rmfield(meta,'f');
% end
%
% % get strip areas and alignment stats for quality selection
% meta.A = cellfun(@(x,y) polyarea(x,y), meta.x,meta.y);
% if isfield(meta,'avg_rmse')
%     meta.scene_alignment_meanrmse= meta.avg_rmse;
% else
%     meta.scene_alignment_meanrmse= cellfun( @(x) nanmean(x.rmse),...
%         meta.scene_alignment);
% end
%
% % find index of tile in tile def database
% tileInd = find(strcmp(tileDefs.I,tileName));
%
% % get tile boundaries with buffer
% x0=tileDefs.x0(tileInd)-buffer;
% y0=tileDefs.y0(tileInd)-buffer;
% x1=tileDefs.x1(tileInd)+buffer;
% y1=tileDefs.y1(tileInd)+buffer;
%
% fprintf('Building water mask tile from %s\n',waterTileDir)
% landTile = getTileWaterMask(waterTileDir,tileName,x0,x1,y0,y1,res);
%
% % make array of subtile boundary coordinates
% subx0=x0:subtileSize:tileDefs.x1(tileInd)-subtileSize-buffer;
% suby0=y0:subtileSize:tileDefs.y1(tileInd)-subtileSize-buffer;
% subx1=tileDefs.x0(tileInd)+subtileSize+buffer:subtileSize:x1;
% suby1=tileDefs.y0(tileInd)+subtileSize+buffer:subtileSize:y1;
%
% [subx0,suby0] = meshgrid(subx0,suby0);
% [subx1,suby1] = meshgrid(subx1,suby1);
%
% subN=numel(subx0);

outNames=dir([outDir,'/',tileName,'_*_10m.mat']);
outNames=cellfun(@(x) [outDir,'/',x], {outNames.name},'uniformoutput',0);


for n=1:length(outNames)
    
    clear x y z offsets za fileNames fileNames0 fa dX dY dZ land offsetErrMax
    
    fprintf('subtile %d of %d\n',n,length(outNames))
    
    outName = outNames{n};
    
    load(outName,'offsetErrMax','dZ')
    
    if ~exist('offsetErrMax','var') || ~exist('dZ','var')
        continue
    end
        
    if ~isinf(offsetErrMax) || ~any(~isnan(dZ) & dZ ~= 0)
        continue
    end
    
    load(outName)
    
    
    % set offsets with horizontal shift failure to zero
    offsets.dx(offsets.dxe == 0) = NaN;
    offsets.dy(offsets.dye == 0) = NaN;
    
    offsets.dxe(isnan(offsets.dx)) = NaN;
    offsets.dye(isnan(offsets.dy)) = NaN;
    
    
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
    
    %number of grid points with coverage
    c0 = any(z,3);
    c0 = sum(c0(:));
    
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
    nsort = nsort(1:i);
    
    if ~any(nsort ~= sort(nsort))
        continue
    end
    
    fprintf('redoing adjustment\n')
    
    nsort=sort(nsort);
    
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
    
    
    fprintf('saving adjustment variables to %s\n',outName)
    save(outName,'dZ','dX','dY','-append');
    
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
    
    % apply a pixel-by-pixel filter to remove outliers
    fa = pairwiseDifferenceFilter(za,'mask',land,'minmad',2);
    fprintf('saving fa to %s\n',outName)
    save(outName,'fa','-append');
    
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

        segs=segs(:);
        r = [r(:);segs(2:end)];
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

