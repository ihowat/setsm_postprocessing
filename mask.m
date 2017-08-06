function  m = mask(varargin)
% MASK ArcticDEM masking algorithm
%
% m = mask(demFile,effectiveBandwidth,abScaleFactor,meanSunElevation) 
% returns the mask stucture (m.x,m.y,m.z) for the demFile using the
% given image parameters.
%
% m = mask(...,maxDigitalNumber,previewPlot) maxDigitalNumber is optional 
% for rescaling the orthoimage to the original source image range. 
% If it's mot included or is empty, no rescaling will be applied. If 
% previewPlot == 'true', a *_maskPreview.tif image will be saved to the
% same directory as the demFile that shows the DEM hillshade with and
% without the mask applied.
%
% m = mask(demFile,meta) returns the mask stucture (m.x,m.y,m.z) for the
% demFile and meta structure, where meta is the output of readSceneMeta.
% Required fields in the meta structure are:
% 'image_1_wv_correct'
% 'image_1_effbw' 
% 'image_1_abscalfact'
% 'image_1_mean_sun_elevation'
% additionally, if image_1_wv_correct==1, the image_1_max field is also
% required.
%
% REQUIRED FUNCTIONS: readGeotiff, DataDensityMap, rescaleDN, edgeSlopeMask
% cloudMask, DG_DN2RAD, waterMask
%
% Ian Howat, ihowat@gmail.com
% 25-Jul-2017 12:49:25

%% Parse inputs
maxDN = [];
previewFlag = false;
if nargin == 2
    demFile=varargin{1};
    meta=varargin{2};
    wv_correctflag=meta.image_1_wv_correct;
    effbw = meta.image_1_effbw;
    abscalfact=meta.image_1_abscalfact;
    mean_sun_elevation = meta.image_1_mean_sun_elevation;
    if wv_correctflag
        maxDN = meta.image_1_max;
    end
elseif nargin >= 4 || nargin <= 6
    demFile=varargin{1};
    effbw =varargin{2};
    abscalfact =varargin{3};
    mean_sun_elevation =varargin{4};
    if nargin == 5 || nargin == 6; maxDN = varargin{5}; end;
    if nargin == 6; previewFlag = varargin{6}; end;
else
    error('2, 4, 5 or 6 input arguments required')
end

% intialize output
m = [];

% Ortho and Matchtag Files
orthoFile = strrep(demFile,'dem.tif','ortho.tif');
mtFile= strrep(demFile,'dem.tif','matchtag.tif');

%% read data and make initial arrays

% read DEM data
z = readGeotiff(demFile);

% initialize output
m.x = z.x;
m.y = z.y;
m.z = false(size(z.z));
m.info = z.info;
m.Tinfo = z.Tinfo;

% parse structure
x = z.x;
y = z.y;
z = z.z;
z(z == -9999) = NaN;

% orignal size for rescaling
sz0=size(z);

% original background for rescaling
bg0 = isnan(z);

% downscale to 8m
z = imresize(z,.25);
x = imresize(x,.25);
y = imresize(y,.25);

% read matchtag
mt = readGeotiff(mtFile);
mt = mt.z;

% make data density map
P = DataDensityMap(mt,21);

% downscale to 8m
P = imresize(P,0.25);

% set P no data
P(isnan(z)) = NaN;

%read ortho
or = readGeotiff(orthoFile);
or = or.z;

% data re-scaling
if ~isempty(maxDN)
    fprintf('rescaled to: 0 to %d\n',maxDN)
    or=rescaleDN(or,maxDN);
end

%convert to radiance
[~,satID]=fileparts(demFile);
satID=upper(satID(1:4));
or = DG_DN2RAD(or, satID,effbw, abscalfact);

fprintf('radiance value range: %.2f to %.2f\n',...
    min(or(:)),max(or(:)))

%rescale ortho to 8m
or = imresize(or,.25);

%% edge crop
M = edgeSlopeMask(x,y,z);

% apply mask
z(~M) = NaN;

% data existence check
if ~any(~isnan(z(:))); return; end

clear M

%% Water Mask
% set no data z values to 0 in ortho
or(isnan(z)) = 0;
P(isnan(z)) = 0;
M = waterMask(or,mean_sun_elevation,P);

% apply mask
z(~M) = NaN;
P(~M) = 0;

% data existence check
if ~any(~isnan(z(:))); return; end

%% Cloud Filter
M = cloudMask(z,or,P);

% apply mask
z(M) = NaN;

%% finalize mask
M = ~isnan(z);

% data existence check
if ~any(M(:)); return; end

M = bwareaopen(M,500);

M = imresize(M,sz0,'nearest');
M(bg0) = false;

% add to structure
m.z = M;
% 

if ~previewFlag; return; end

%% Display output
 % 
 % read DEM data
z = readGeotiff(demFile);
x = z.x;
y = z.y;
z = z.z;
z(z == -9999) = NaN;
 
z_masked = z;
z_masked(~M) = NaN;
 
rf = 0.25;
z = imresize(z,rf);
x = imresize(x,rf); % only needed for hillshade
y = imresize(y,rf); % only needed for hillshade
z_masked = imresize(z_masked,rf);
 
hill = hillshade(z,x,y);
hill_masked = hillshade(z_masked,x,y);
 
h=figure(1);
set(gcf,'color','w','visible','off')
hpos = get(gcf,'position');
hpos(3) = 2000;
hpos(4) = 1000;
set(gcf,'position',hpos);
subplot(1,2,1)
imagesc(hill,'alphadata',~isnan(z)); colormap gray; axis equal
subplot(1,2,2)
imagesc(hill_masked,'alphadata',~isnan(z_masked)); colormap gray; axis equal
print(h,'-dtiff','-r100',strrep(demFile,'_dem.tif','_maskPreview.tif'));
close(h)
