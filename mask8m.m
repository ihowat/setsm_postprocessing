function  m = mask8m(varargin)
% MASK masking algorithm for 8m resolution data
%
% m = mask(demFile) 
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
%
% REQUIRED FUNCTIONS: readGeotiff, DataDensityMap,  edgeSlopeMask,
% DataDensityMask
%
% Ian Howat, ihowat@gmail.com
% 25-Jul-2017 12:49:25

%% Parse inputs
previewFlag = false;
demFile=varargin{1};

% intialize output
m = [];
mtFile= strrep(demFile,'dem.tif','matchtag.tif');

%% read data and make initial arrays

% read DEM data
z = readGeotiff(demFile);

% read matchtag
mt = readGeotiff(mtFile);

% size consistency checks
if any(size(z.z) ~= size(mt.z))
    error('size of dem and matchtag rasters dont match')
end

% initialize output
m.x = z.x;
m.y = z.y;
m.z = false(size(z.z));
m.info = z.info;
m.Tinfo = z.Tinfo;

% parse structures
x = z.x;
y = z.y;
z = z.z;
z(z == -9999) = NaN;
mt = mt.z;

% make data density map
P = DataDensityMap(mt,21);

% set P no data
P(isnan(z)) = NaN;

%% edge crop
M = edgeSlopeMask(x,y,z);

M(isnan(z)) = false;

% data existence check
if ~any(M(:)); return; end

clear z

%% Data Density Filter
n= 21;
Pmin=0.90;
Amin=1000;
Amax=1000;

M =DataDensityMask(M,'n',n,'Pmin',Pmin,'Amin',Amin,'Amax',Amax);

% data existence check
if ~any(M(:)); return; end

M = bwareaopen(M,500);

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
