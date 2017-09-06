function  m = remaMask2a(demFile)
% remaMask2a mask dem using point density and slopes deviations
%
% m = remaMask2a(demFile) masks the dem file by first applying an edgemask
% and then using an iterative bad pixel search starting from points of low
% point density and high standard deviation in slope.

% plot hillshade browse image?
plotBrowse = false;

%% Filenames
mtFile= strrep(demFile,'dem.tif','matchtag.tif');

%% read DEM data
z = readGeotiff(demFile);

%% initialize output
m.x = z.x;
m.y = z.y;
m.z = false(size(z.z));
m.info = z.info;
m.Tinfo = z.Tinfo;

x = z.x;
y = z.y;
z = z.z;
z(z == -9999) = NaN;

%% edge crop

% slope of z
[dx,dy] = gradient(z,x,y);

% aspect and grade (aspect is used in the cloud filter)
[~,r] = cart2pol(dx,dy);

% meang grade over n-pixel kernel
n=5;
Mr = conv2(r,ones(n)/(n.^2),'same');

% mask mean slopes greater than 1
M = Mr < 1;

% dilate high mean slope pixels by 51 pixels and set to false
M(imdilate(Mr > 1,ones(13))) = false;

clear Mr

% fill interior holes since we're just looking for edges here.
M = imfill(M,'holes');

% get rid of isolated little clusters of data
M = bwareaopen(M, 1000);

% no data check
if ~any(M(:)); fprintf('boundary filter removed all data\n'); return; end;

% find data coverage boundaries
B = bwboundaries(M, 8, 'noholes');

% merge all the data clusters
B = cell2mat(B);

% find outer data boundary
k = boundary(B(:,2),B(:,1),0.5); %allows for some inward "bending" to conform to the data.
%k = convhull(B(:,2),B(:,1));% straight line trace of the outer edge.

% build edge mask
M = poly2mask(B(k,2),B(k,1), size(M,1),size(M,2));

clear B k

% apply mask
z(~M) = NaN;

clear M

dx(isnan(z))=NaN;
dy(isnan(z))=NaN;


%% Iterative expanding mt density/slope mask

% standard deviation of slopes
n=5;
Mdx =conv2(dx,ones(n)/(n.^2),'same');
Sdx=sqrt( conv2(dx.^2,ones(n)/(n.^2),'same') - Mdx.^2 );

Mdy =conv2(dy,ones(n)/(n.^2),'same');
Sdy=sqrt( conv2(dy.^2,ones(n)/(n.^2),'same') - Mdy.^2 );

Sd = sqrt(Sdx.^2 + Sdy.^2);

% read matchtag
mt = readGeotiff(mtFile);
mt = mt.z;

% make data density map
P = DataDensityMap(mt,n);

% set P no data
P(isnan(z)) = NaN;

%locate probable cloud pixels
M = Sd > 0.75 & P < 1;

% remove small clusters
M = bwareaopen(M,10);

% initialize masked pixel counters
N0 = sum(M(:));
N1 = inf;

% background mask
MM = isnan(z) | isnan(Sd);

% expand mask to surrounding bad pixels, stop when mask stops growing
while N0 ~= N1
    
    N0 = N1; % set new to old
    M = imdilate(M,ones(11)); % dialate the mask
    M(MM | Sd < 0.1) = false; % unmask low sd pixels
    N1 = sum(M(:)); % count the number of new masked pixels.
    
end

% remove small data gaps
M = ~bwareaopen(~M,5000);

% remove border effect
M = M | imdilate(isnan(z),ones(5));

% remove small data gaps
M = ~bwareaopen(~M,5000);

% put mask into output stucture with good data == 1
m.z = ~M;


if ~plotBrowse; return; end


%% Plot browse hillshade image

z = readGeotiff(demFile);
x = z.x;
y = z.y;
z = z.z;
z(z == -9999) = NaN;
hill = hillshade(z,x,y);


h=figure(1);
%set(gcf,'color','w','visible','off')
set(gcf,'color','w','visible','on')
hpos = get(gcf,'position');
hpos(3) = 2000;
hpos(4) = 1000;
set(gcf,'position',hpos);
subplot(1,2,1)
imagesc(hill,'alphadata',~isnan(z)); colormap gray; axis equal
subplot(1,2,2)
imagesc(hill,'alphadata',m.z); colormap gray; axis equal
%print(h,'-dtiff','-r100',strrep(demFile,'_dem.tif','_mask3e_8m_preview.tif'));
input('')
close(h)


