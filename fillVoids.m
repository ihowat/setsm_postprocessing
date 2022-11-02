function fillVoids(subTileName)
%fillVoids
% subTileName example: 49_10_1_1

sl = split(subTileName,'_');
if startsWith(sl{1},'utm')
    tileName = strjoin(sl(1:3),'_');
else
    tileName = strjoin(sl(1:2),'_');
end
projstr = 'polar stereo north';
demDir = 'V:\pgc\data\elev\dem\setsm\ArcticDEM\mosaic\v4.1\2m';
demFile = dir([demDir,'/',tileName,'/',subTileName,'_2m_dem.tif']);
demFile = cellfun(@(x) [demDir,'/',tileName,'/',x], {demFile.name}, 'uniformOutput', false);
[tileProjName,projstr] = getProjName(tileName,projstr);

fprintf('Loading tile %s water mask\n',tileName)
waterTileDir='V:\pgc\data\projects\arcticdem\watermasks\global_surface_water\tiled_watermasks';
waterMask = dir([waterTileDir,'/',tileName,'_water.tif']);
waterMask = cellfun(@(x) [waterTileDir,'/',x], {waterMask.name}, 'uniformOutput',false);
if isempty(waterMask)
    fprintf('Tile %s water mask file does not exist in %s\n',tileName,waterTileDir)
    return
end
waterMaskFile = waterMask{1};

% load strip database
databaseFile = 'V:\pgc\data\scratch\claire\repos\setsm_postprocessing_pgc\arcticDEMdatabase4_2m_v4_20200806.mat';
changePath= 'V:/pgc'; %if set, will change the path to the directory from what's in the database file. set to [] if none.

meta=load(databaseFile,'fileName','x','y');

voidMaskFile=strrep(demFile,'.tif','_voidMask.mat');
outName=strrep(demFile,'_dem.tif','_voidFilled_dem.tif');

%% load data
fprintf('loading files\n')

if ~exist(voidMaskFile,'file')
    error('%s doesnt exist\n',voidMaskFile)
end

dem = readGeotiff(demFile);
N= readGeotiff(strrep(demFile,'_dem','_N'));
dem.N= N.z; clear N;
z_mad=readGeotiff(strrep(demFile,'_dem','_mad'));
dem.z_mad=z_mad.z; clear z_mad;
tmax=readGeotiff(strrep(demFile,'_dem','_tmax'));
dem.tmax=tmax.z; clear tmax
tmin=readGeotiff(strrep(demFile,'_dem','_tmin'));
dem.tmin=tmin.z; clear tmin;

land = readGeotiff(waterMaskFile);
land.z = land.z== 0;
land = interp2(land.x,land.y(:),single(land.z),dem.x,dem.y(:),'*nearest');
land(isnan(land)) = 0;
land = logical(land);

%load(coastlinePolyFile);
load(voidMaskFile)

% convert nodata to nans
dem.z(dem.z == -9999) = NaN;

% apply voidMask to dem.z
dem.z(~voidMask) = NaN;

%% build land mask
% make polyshape of this tile with buffer to ensure coverage of border
% cells
% res=100;
% x0  = min(dem.x);
% x1  = max(dem.x);
% y0  = min(dem.y);
% y1  = max(dem.y);
% tilePoly = polyshape([x0-res;x0-res;x1+res;x1+res],[y0-res;y1+res;y1+res;y0-res]);
%
% land =  false(length(dem10.y),length(dem10.x));
%
% for i=1:length(coastlinePoly)
%     if overlaps(tilePoly,coastlinePoly(i))
%
%         landPoly = intersect(tilePoly,coastlinePoly(i));
%
%         NR = [0;find(isnan(landPoly.Vertices(:,1)));...
%             length(landPoly.Vertices(:,1))+1];
%
%         for nr = 1:length(NR)-1
%             xp = landPoly.Vertices(NR(nr)+1:NR(nr+1)-1,1);
%             yp = landPoly.Vertices(NR(nr)+1:NR(nr+1)-1,2);
%             land(roipoly(dem.x,dem.y,land,xp,yp))=true;
%         end
%
%     end
% end

%% loop through void polys to select strips to use for fill
f1=figure;
i=1;
for i=1:length(voidPolys)

    fprintf('void %d of %d\n',i,length(voidPolys))

    % check to make sure void area is more than a few pixels
    if area(voidPolys(i)) < 500
        continue
    end

    % indeices of stirps overlaping this void
    ind=stripSearch(meta.x,meta.y,voidPolys(i));

    if isempty(ind)
        fprintf('no strips overlap, skipping\n')
        continue
    end

    % convert meta names into dem names
    fileNames =  strrep(meta.fileName(ind),'meta.txt','dem.tif');
    browseNames = strrep(meta.fileName(ind),'meta.txt','dem_10m_shade.tif');

    if ismac
        fileNames  = strrep(fileNames,'/fs/project/howat.4','/Users/ihowat/project');
        browseNames  = strrep(browseNames,'/fs/project/howat.4','/Users/ihowat/project');
    else
        fileNames  = strrep(fileNames,'/mnt/pgc','V:/pgc');
        browseNames  = strrep(browseNames,'/mnt/pgc','V:/pgc');
    end

    % rectangular coordinate range of polygon with buffer
    buffer = 100;
    x0=min(voidPolys(i).Vertices(:,1)) - buffer;
    x1=max(voidPolys(i).Vertices(:,1))  + buffer;
    y0=min(voidPolys(i).Vertices(:,2))  - buffer;
    y1=max(voidPolys(i).Vertices(:,2))  + buffer;

    % if region too small, make larger to see more clearly
    if (x1 - x0) < 2000
        xm = (x0 + x1)./2;
        x0 = xm - 1000;
        x1 = xm + 1000;
    end
    if (y1 - y0) < 2000
        ym = (y0 + y1)./2;
        y0 = ym - 1000;
        y1 = ym + 1000;
    end

    % initialize strip selection vector, 1 means need to decide
    usevec= ones(size(browseNames));

    % loop through strips to examine and select or discard
    % keep looping until all decided on (no 1's in usevec)
    j=1;
    while any(usevec == 1)

        % if cycled through, start at beginning
        if j > length(usevec)
            j = 1;
        end

        % skip strips already rejected or selected
        if usevec(j) == 0 || usevec(j) == 2
            j=j+1;
            continue
        end

        fprintf('dem strip %d of %d\n',j,length(browseNames))
        z1 = readGeotiff(browseNames{j},'map_subset',[x0 x1 y0 y1]);

        if length(z1.x) < 10 || length(z1.y) < 10 || ~any(z1.z(:))
            fprintf('No data,skipping\n')
            usevec(j)=0;
            j=j+1;
            continue
        end

        % plot image with poly
        imagesc(z1.x,z1.y,z1.z)
        axis xy equal; colormap gray
        hold on
        plot(voidPolys(i))

        %get input
        stdin = input('(u)se, (r)eject, (h)old, (s)top?\n','s');

        switch lower(stdin)
            case 'r' % reject this strip
                usevec(j) = 0;
            case 'h'  % hold in the deck
                usevec(j) = 1;
            case 'u' % use this strip
                usevec(j) = 2;
            case 's' % stop selecting
                usevec(usevec ~= 2) = 0;
                clf

                break
            otherwise % unrecognized input, do over
                clf
                continue
        end

        clf

        j=j+1;

    end

    % turn usevec into index of selected dems
    usevec = usevec==2;

    if ~any(usevec)
        fprintf('no strips accepted, leaving void\n')
        continue
    end

    % crop filename vector to those selected
    fileNames=fileNames(usevec);

    % make date vector
    [~,name] =  cellfun(@fileparts,fileNames,'uniformoutput',0);
    t=cellfun(@(x) datenum(parsePairnameDatestring(x),'yyyymmdd'),name)';

    % re-specify rectangular sampling area
    buffer = 100;
    x0=min(voidPolys(i).Vertices(:,1)) - buffer;
    x1=max(voidPolys(i).Vertices(:,1))  + buffer;
    y0=min(voidPolys(i).Vertices(:,2))  - buffer;
    y1=max(voidPolys(i).Vertices(:,2))  + buffer;

    fprintf('extracting %d strip subsets\n',length(fileNames))
    [x,y,z] =extractSubGrid(fileNames,x0,x1,y0,y1,2);

    % make sure more than one dem for coregistration/adustment
    if size(z,3) > 1

        %make a vector of z's that ar belonging to the same strip
        [~,stripid] =  cellfun(@fileparts,fileNames,'uniformoutput',0);
        stripid =  cellfun(@(x) x(1:47),stripid,'uniformoutput',0);
        [unique_stripids,~,strip_ind] = unique(stripid);

        % make sure more than one unique strip
        if length(unique_stripids) > 1

            landsub = interp2(dem.x,dem.y(:),single(land),x,y,'*nearest');
            landsub(isnan(landsub)) = 0;
            landsub = logical(landsub);

            fprintf('performing pairwise coregistration, ')
            offsets=coregisterStack(x,y,z,landsub,strip_ind);

            offsets.dx(offsets.dxe == 0) = NaN;
            offsets.dy(offsets.dye == 0) = NaN;

            offsets.dxe(isnan(offsets.dx)) = NaN;
            offsets.dye(isnan(offsets.dy)) = NaN;

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
            z= za;
            clear za

            fa = pairwiseDifferenceFilter(z,'mask',landsub,'minmad',2);

            z(~fa) = NaN;
        end

    end

    N = uint8(sum(~isnan(z),3));
    z_mad = mad(z,1,3);

    t=t-datenum('1/1/2000 00:00:00');
    t=reshape(t,1,1,[]);
    t=repmat(t,size(z));
    t(isnan(z))=NaN;
    tmax = max(t,[],3);
    tmin = min(t,[],3);
    %tmean = nanmean(t,3);

    tmax = uint16(tmax);
    tmin = uint16(tmin);
    z = nanmedian(z,3);


    % make sure x and y are round for locating in master grid
    x=round(x);
    y=round(y);

    % get corner indexes for this grid in the master
    col0  = find(x(1) == dem.x);
    if isempty(col0)
        col0 = 1;
        crop=find(dem.x(1) == x);
        x = x(crop:end);
        z = z(:,crop:end);
    end
    col1  = find(x(end) == dem.x);
    if isempty(col1)
        col1 = length(dem.x);
        crop=find(dem.x(end) == x);
        x = x(1:crop);
        z = z(:,1:crop);
    end
    row0  = find(y(1) == dem.y);
    if isempty(row0)
        row0 = 1;
        crop=find(dem.y(1) == y);
        y = y(crop:end);
        z = z(crop:end,:);
    end
    row1  = find(y(end) == dem.y);
    if isempty(row1)
        row1  = length(dem.y);
        crop=find(dem.y(end) == y);
        y = y(1:crop);
        y = y(1:crop,:);
    end

    % subset the current mosaic by range of subtile
    z1 = dem.z(row0:row1,col0:col1);
    land1 = land(row0:row1,col0:col1);
    N1 = dem.N(row0:row1,col0:col1);
    z_mad1 = dem.z_mad(row0:row1,col0:col1);
    tmax1 = dem.tmax(row0:row1,col0:col1);
    tmin1 = dem.tmin(row0:row1,col0:col1);

    % find overlapping non-nan and non-water pixels
    n_overlap = ~isnan(z1(:)) & ~isnan(z(:)) & land1(:);

    % check if overlapping pixels exist to determine if blending is needed
    if any(n_overlap)

        %get median difference between this subtile and mosaic subset
        dz_med = median(z(n_overlap) - z1(n_overlap));

        % subtile number-weighted alignement
        z = z - dz_med; % shift the fill dem to the mosaic

        % fill missing data in this subtile with current data in mosaic
        % subset
        z(isnan(z) & ~isnan(z1)) =...
            z1(isnan(z) & ~isnan(z1));

        z_mad(isnan(z_mad) & ~isnan(z_mad1)) =...
            z_mad1(isnan(z_mad) & ~isnan(z_mad1));


        % Create the blending array by setting zeros at the far edge of
        % subtile/tile overlap and ones at the other edge, linearly
        % interpolating between the two.
        % find where data is missing in both subtile and tile
        buffA = single(~(~isnan(z) & ~isnan(z1)));

        % set pixels with data in both subtile and tiles to NaN as an
        % interpolation flag for inpaint_nans
        buffA(~buffA) = NaN;

        % set boundaries of blend array to zero
        buffA(1,isnan(buffA(1,:))) = 0;
        buffA(end,isnan(buffA(end,:))) = 0;

        buffA(isnan(buffA(:,1)),1) = 0;
        buffA(isnan(buffA(:,end)),end) = 0;

        % interpolate linearly across NaNs
        buffA=inpaint_nans(double(buffA),2);

        % find where there is data in the tile subset
        notMissing = ~isnan(z1);

        % blend the subtile and mosaic subset, applying the edge-distance
        % weighting
        z(notMissing) = z(notMissing).*buffA(notMissing) +...
            z1(notMissing).*(1- buffA(notMissing));


        z_mad(notMissing) = z_mad(notMissing).*buffA(notMissing) +...
            z_mad1(notMissing).*(1- buffA(notMissing));


    else

        % if no pixels overlap, just add the subset data into the subtile,
        % replacing the NaNs
        z(~isnan(z1(:))) = z1(~isnan(z1(:)));

        z_mad(~isnan(z1(:))) = z_mad1(~isnan(z1(:)));

    end

    % place the belended substile into the tile
    dem.z(row0:row1,col0:col1) = z;
    dem.z_mad(row0:row1,col0:col1) = z_mad;

    N1(~notMissing) = N(~notMissing);
    dem.N(row0:row1,col0:col1) = N1;

    tmax1(~notMissing) = tmax(~notMissing);
    dem.tmax(row0:row1,col0:col1) = tmax1;

    tmin1(~notMissing) = tmin(~notMissing);
    dem.tmin(row0:row1,col0:col1) = tmin1;

end


dem.z(isnan(z)) = -9999;

writeGeotiff(outName,x,y,dem.z,4,-9999,projstr)

outNameTif = strrep(outName,'_dem.tif','_N.tif');
writeGeotiff(outNameTif,x,y,dem.N,1,0,projstr)

z_mad(isnan(z_mad)) = -9999;
outNameTif = strrep(outName,'_dem.tif','_mad.tif');
writeGeotiff(outNameTif,x,y,dem.z_mad,4,-9999,projstr)

outNameTif = strrep(outName,'_dem.tif','_tmax.tif');
writeGeotiff(outNameTif,x,y,dem.tmax,2,0,projstr)

outNameTif = strrep(outName,'_dem.tif','_tmin.tif');
writeGeotiff(outNameTif,x,y,dem.tmin,2,0,projstr)
