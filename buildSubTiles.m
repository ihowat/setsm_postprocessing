function buildSubTiles(tileName)

%tileName='47_13';
res=10; % ouput mosaic resolution in meters
subtileSize=1000; % subtile dimensions in meters
buffer=100; % size of tile/subtile boundary buffer in meters
redoFlag = 4; % flag for treating existing subtile .mat files: 1=skip, 2=redo coregistration, 3=redo adjustment
maxNumberOfStrips=100; % maximum number of strips to load in subtile

%paths/files
if ismac
    tileDefFile = 'PGC_Imagery_Mosaic_Tiles_Arctic.mat';
    databaseFile = '/Users/ihowat/gdrive/projects/earthdem/earthdem_database_unf.mat';
    outDir = ['/Users/ihowat/project/howat.4/earthdem/earthdem_mosaic_testing_1km/',tileName];
    addpath('/Users/ihowat/unity-home/demtools');
    coastlinePolyFile='gshhg_237_alaska_coastline_3413.mat';
    lakePolyFile='gshhg_237_alaska_lakes_3413.mat';
else
    tileDefFile = 'PGC_Imagery_Mosaic_Tiles_Arctic.mat'; %PGC/NGA Tile definition file
     databaseFile = 'arcticdem_database_unf_pgcpaths.mat';
     outDir = ['/mnt/pgc/data/scratch/claire/pgc/arcticdem/mosaic/2m_v4/',tileName,'/subtiles'];
     %addpath('/home/howat.4/demtools');
     waterTileDir='/mnt/pgc/data/scratch/claire/pgc/arcticdem/coastline/water_tiles';
     %coastlinePolyFile='gshhg_237_alaska_coastline_3413.mat';
     %lakePolyFile='gshhg_237_alaska_lakes_3413.mat';

end

if ~exist(outDir,'dir')
    mkdir(outDir)
end

fprintf('Loading tile definition file\n')
tileDefs=load(tileDefFile);
fprintf('Loading strip database and getting tile overlaps\n')
meta=load(databaseFile);

meta.A = cellfun(@(x,y) polyarea(x,y), meta.x,meta.y);
meta.scene_alignment_meanrmse= cellfun( @(x) mean(x.rmse), meta.scene_alignment);

tileInd = find(strcmp(tileDefs.I,tileName));

% get tile boundaries with buffer
x0=tileDefs.x0(tileInd)-buffer;
y0=tileDefs.y0(tileInd)-buffer;
x1=tileDefs.x1(tileInd)+buffer;
y1=tileDefs.y1(tileInd)+buffer;

% load coastline tile polyshape
fprintf('Loading tile %s coastline\n',tileName)
coastlinePolyFile = dir([waterTileDir,'/',tileName,'_coast.mat']);
coastlinePolyFile = cellfun(@(x) [waterTileDir,'/',x], {coastlinePolyFile.name}, 'uniformOutput',false);
if isempty(coastlinePolyFile)
    fprintf('Tile %s coastline file does not exist in %s\n',tileName,waterTileDir)
    return
end
load(coastlinePolyFile{1});

% load lakes tile polyshape
fprintf('Loading tile %s lakes\n',tileName)
lakePolyFile = dir([waterTileDir,'/',tileName,'_lakes.mat']);
lakePolyFile = cellfun(@(x) [waterTileDir,'/',x], {lakePolyFile.name}, 'uniformOutput',false);
%lakePolyFile{1}='gshhg_237_alaska_lakes_3413.mat';
if isempty(lakePolyFile)
    fprintf('Tile %s lakes file does not exist in %s\n',tileName,waterTileDir)
    waterPoly(1)=polyshape();
else
    load(lakePolyFile{1});
    waterPoly = lakePoly;
end

% make polyshape of this tile with buffer to ensure coverage of border
% cells
tilePoly = polyshape([x0-res;x0-res;x1+res;x1+res],[y0-res;y1+res;y1+res;y0-res]);

% make polygons of land and water over this tile
i=1;
count=1;
clear landPoly
for i=1:length(coastlinePoly)
    if overlaps(tilePoly,coastlinePoly(i))
        landPoly(count) = intersect(tilePoly,coastlinePoly(i));
        count=count+1;
    end
end

if isempty(landPoly)
    fprintf('No land surface in tile %s, returning\n',tileName)
end
landPoly = union(landPoly);
clear coastlinePoly

% % polygons of lakes over this tile
% i=1;
% count=1;
% clear waterPoly
% waterPoly(1)=polyshape();
% for i=1:length(lakePoly)
%     if overlaps(tilePoly,lakePoly(i))
%         waterPoly(count) = intersect(tilePoly,lakePoly(i));
%         count=count+1;
%     end
% end

if ~isempty(waterPoly)
    waterPoly = union(waterPoly);
end
clear lakePoly

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
    
    outName = [outDir,'/',tileName,'_',num2str(n),'_10m.mat'];
    
    if exist(outName,'file')
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
                
        end
    end
    
    % find overlapping strips
    if ~exist('fileNames','var')
        
        % make sure this subtile has land surface in it
        subtilePoly = polyshape([subx0(n);subx0(n);subx1(n);subx1(n)],...
            [suby0(n);suby1(n);suby1(n);suby0(n)]);
        
        if ~overlaps(subtilePoly,landPoly)
            fprintf('No land in subtile %d, skipping\n',n)
            continue
        end
        
        % make land surface polygon within this subtile for searching
        sublandPoly = intersect(subtilePoly,landPoly);
        
        ind=stripSearch(meta.x,meta.y,sublandPoly);
        % ind=stripSearch(meta.x,meta.y,subx0(n),subx1(n),suby0(n),suby1(n));
        
        % check for maximum #'s of overlaps
        if length(ind) > maxNumberOfStrips   
            %first remove single subscenes
            ind(meta.scene_alignment_meanrmse(ind) == 0 |...
                isnan(meta.scene_alignment_meanrmse(ind))) = [];
        end
        
        if length(ind) >  maxNumberOfStrips
            % sort by weigthed area and coverage
            scene_alignment_meanrmse_scaled= -(meta.scene_alignment_meanrmse(ind)-mean(meta.scene_alignment_meanrmse(ind)))./std(meta.scene_alignment_meanrmse(ind));
            A_scaled = (meta.A(ind)-mean(meta.A(ind)))./std(meta.A(ind));
            
            scene_alignment_meanrmse_scaled = scene_alignment_meanrmse_scaled + abs(min(scene_alignment_meanrmse_scaled));
            A_scaled = A_scaled + abs(min(A_scaled));
            
            [~,ind_rank] = sort(A_scaled.*scene_alignment_meanrmse_scaled,'descend');
            
            ind = ind(ind_rank(1:maxNumberOfStrips));
        end
        
        if isempty(ind)
            fprintf('no overlapping strips, skipping\n')
            continue
        end
        
        % convert meta names into dem names
        fileNames  = strrep(meta.fileName(ind),'meta.txt','dem_10m.tif');
        
        if ismac
            fileNames  = strrep(fileNames,'/fs/project/howat.4','/Users/ihowat/project');
        end
        
        % make date vector
        [~,name] =  cellfun(@fileparts,fileNames,'uniformoutput',0);
        t=cellfun(@(x) datenum(x(6:13),'yyyymmdd'),name)';
        
    end
    
    if ~exist('z','var')
        
        fprintf('extracting %d strip subsets ',length(fileNames))
        [x,y,z] =extractSubGrid(fileNames,subx0(n),subx1(n),...
            suby0(n),suby1(n),res);
        
        fprintf('saving x,y,z,t,fileNames to %s\n',outName)
        save(outName,'x','y','z','t','fileNames','-v7.3');
        
    end
    
    
    if ~exist('offsets','var')
        % make land mask
        subtilePoly = polyshape([subx0(n)-res;subx0(n)-res;subx1(n)+res;subx1(n)+res],...
            [suby0(n)-res;suby1(n)+res;suby1(n)+res;suby0(n)-res]);
        
        if ~overlaps(subtilePoly,landPoly)
            fprintf('No land in subtile %d, skipping\n',n)
        end
        
        sublandPoly = intersect(subtilePoly,landPoly);
        
        NR = [0;find(isnan(sublandPoly.Vertices(:,1)));...
            length(sublandPoly.Vertices(:,1))+1];
        
        land =  false(length(y),length(x)); % land mask
        
        for nr = 1:length(NR)-1
            xp = sublandPoly.Vertices(NR(nr)+1:NR(nr+1)-1,1);
            yp = sublandPoly.Vertices(NR(nr)+1:NR(nr+1)-1,2);
            land(roipoly(x,y,land,xp,yp))=true;
        end
        
        if overlaps(subtilePoly,waterPoly)
            subwaterPoly = intersect(subtilePoly,waterPoly);
            NR = [0;find(isnan(subwaterPoly.Vertices(:,1)));...
                length(subwaterPoly.Vertices(:,1))+1];
            for nr = 1:length(NR)-1
                xp = subwaterPoly.Vertices(NR(nr)+1:NR(nr+1)-1,1);
                yp = subwaterPoly.Vertices(NR(nr)+1:NR(nr+1)-1,2);
                land(roipoly(x,y,land,xp,yp))=false;
            end
        end
        
        if ~any(land(:))
            fprintf('No land in subtile %d, skipping\n',n)
        end
        
        fprintf('saving land to %s\n',outName)
        save(outName,'land','-append');
        
        %make a vector of z's that ar belonging to the same strip 
        [~,stripid] =  cellfun(@fileparts,fileNames,'uniformoutput',0);
        stripid =  cellfun(@(x) x(1:47),stripid,'uniformoutput',0);
        [~,~,strip_ind] = unique(stripid);
        
        fprintf('performing pairwise coregistration, ')
        offsets=coregisterStack(x,y,z,land,strip_ind);
        
        fprintf('saving offsets to %s\n',outName)
        save(outName,'offsets','-append');
        
    end

    if ~exist('dZ','var')
        
        if isempty(offsets)
            fprintf('no offsets returned, skipping\n')
            continue
        end
        
        offsets.dx(offsets.dxe == 0) = NaN;
        offsets.dy(offsets.dye == 0) = NaN;
        
        offsets.dxe(isnan(offsets.dx)) = NaN;
        offsets.dye(isnan(offsets.dy)) = NaN;
        
        if ~any(~isnan(offsets.dx))
            fprintf('no offsets returned, skipping\n')
            continue
        end
        
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
            
            [dZ,dX,dY] = adjustOffsets(offsets,'offsetErrMax',offsetErrMax,'min_sigma_dz_coregMax',min_sigma_dz_coregMax,...
                'min_abs_mean_dz_coregMax',min_abs_mean_dz_coregMax,'min_abs_median_dz_coregMax',min_abs_median_dz_coregMax);
            
            %check for coverage of z layers with solutions
            c1 = any(z(:,:,~isnan(dZ)),3);
            c1 = sum(c1(:));
            
            if c1 == c0
                break
            end
            
            it = it+0.1;
        end
        
        
        fprintf('saving adjustment variables to %s\n',outName)
        save(outName,'dZ','dX','dY','offsetErrMax','min_sigma_dz_coregMax',...
            'min_abs_mean_dz_coregMax','min_abs_median_dz_coregMax','-append');
        
    end
    
    if ~exist('za','var')
        
        % layers with missing adjustments
        n_missing = isnan(dZ);
        
        % make adjusted z array
        za=nan(size(z),'single');
        
        % index of layers with adjustments
        iterVec = 1:size(z,3);
        iterVec(n_missing) = [];
        
        fprintf('applying adjustment\n')
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
    make2m(fileNames,t,x,y,dZ,dX,dY,land,fa,outName2m);
         
end

function make2m(fileNames,t,x,y,dZ,dX,dY,land,fa,outName)

% layers with missing adjustments
n_missing = isnan(dZ);

dZ(n_missing) = [];
dX(n_missing) = [];
dY(n_missing) = [];
fileNames(n_missing) = [];
fa(:,:,n_missing) = [];
t(n_missing) = [];

fileNames = strrep(fileNames,'_10m.tif','.tif');

[x,y,z] =extractSubGrid(fileNames,min(x),max(x),...
    min(y),max(y),2);

% make adjusted z array
za=nan(size(z),'single');

for k=1:size(z,3)
    zak = interp2(x + dX(k),y + dY(k), z(:,:,k) + dZ(k),...
        x,y,'*linear');
    
    if any(any(~fa(:,:,k)))
        
        zak(~imresize(fa(:,:,k),size(za(:,:,k)),'nearest'))=NaN;
    end
    
    za(:,:,k) = zak;
    
end

za_med = nanmedian(za,3);
%za_std =  nanstd(za,[],3);
za_mad = mad(za,1,3);
N = uint8(sum(~isnan(za),3));

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
save(outName,'x','y','za_med','land','za_mad','N','tmax','tmin');



