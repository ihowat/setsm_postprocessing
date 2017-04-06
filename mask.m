function m = mask(demFile)
%m = mask(demFile) create ArcticDEM scene mask
%
% m = mask(demFile) returns the data mask structure m.x, m.y, m.z for the
%dem in demFile. Mask uses the DEM, orthoimage and matchtag files, which
%are assumed to be int he same path whith the same basename and standard
%suffixes. Mask combines match point-density and orthoimage entropy to
%locate large water bodies and uses the match point-density and the
%standard deviation of slope to locate both water and cloads
%
% Ian Howat, ihowat@gmail.com
% Version 1. 06-Apr-2017 14:14:36

% make file names
mtFile = strrep(demFile,'dem.tif','matchtag.tif');
orFile = strrep(demFile,'dem.tif','ortho.tif');
metaFile = strrep(demFile,'dem.tif','meta.txt');


% There is some sort of scaling difference between the wv_correct and non
% wv_correct - applied othos that screws up the conversion to uint8 in
% entropyfilt (which uses im2uint8) unless I convert to uint8 first. So
% need to check if wvc applied.
wvcFlag=false;

% check if meta file exists
if ~exist(metaFile,'file')
   warning('no meta file found, assuming no wv_correct applied\n');
else
    % exists, so read meta file
    c=textread(metaFile,'%s','delimiter','\n');
    
    %find SETSM version
    str='Image 1 wv_correct=';
    r=find(~cellfun(@isempty,strfind(c,str)));
    value=deblank(strrep(c{r(1)},str,''));
    wvcFlag = logical(str2num(value));
    
    if wvcFlag
        fprintf('wv_correct applied\n')
    else
        fprintf('wv_correct not applied\n');
    end
end

%% Data density image

% read the matchtag and ortho files
mt=readGeotiff(mtFile);

% calculate data density map using default (11x11 kernel)
P = DataDensityMap(mt.z);

% resize to 8m - resizing the mt directly results in point loss
P = imresize(P,0.25);

%% Entropy Image

% read ortho
or=readGeotiff(orFile);
or = or.z; % get rid of map data

% resize ortho to 8m
or = imresize(or,0.25);

% subtraction image
or_subtraction =  movmax(or,5) - movmin(or,5);
 
or_subtraction(imdilate(or == 0,ones(9))) = 0;

if ~wvcFlag; or_subtraction = uint8(or_subtraction); end

 % entropy image
J = entropyfilt(or_subtraction,ones(5));


%% Entropy & Density Mask

% use a two-minimum threshold, the first is high density, low entropy (e.g.
% snow and ice), the second is low density, higher entropy (e.g. forest
% canopy)
M = (P > 0.9 & J > 0.2) | (P > 0.3 & J > 1) ; 

% just concerned with coastline here so fill the interior data gaps
M = imfill(M,'holes');

% get rid of isolated little clusters of data. we could get rid of more if
% a minumum "island" size is known
M = bwareaopen(M, 1000); 

% expand to ensure boundary/coast coverage
M = imdilate(M, ones(7)); 

% set background to false
M(or ==0) = false;


%% Sigma Slope & Density Mask

% read dem
z = readGeotiff(demFile);

% x = z.x; only needed for hillshade 
% y = z.y; only needed for hillshade
z = z.z;

z(z == -9999) = NaN;

z = imresize(z,.25);
% x = imresize(x,.25); only needed for hillshade
% y = imresize(y,.25); only needed for hillshade

% construct hillshade for debugging
%hill = hillshade(z,x,y);

[sx,sy] = gradient(z,8); % calculate slopes

[~,rho] = cart2pol(sx,sy); % vector slope

k=11; % convolution kernel size
rhomn =conv2(rho,ones(k)/(k.^2),'same'); % mean slope
rhosd=sqrt( conv2(rho.^2,ones(k)/(k.^2),'same') - rhomn.^2 ); % std dev

% sigma slope threshold function
rhosdthresh= (1.*P).^4;

% apply theshold
M1 = rhosd <= rhosdthresh;

% remove data bits < 2000
M1 = bwareaopen(M1,2000);

% remove data gaps < 2000
M1 = ~bwareaopen(~M1,2000);

% build output structure
m.z = imresize(M & M1,size(mt.z),'nearest');
m.x = mt.x;
m.y = mt.y;
m.info = mt.info;
m.Tinfo = mt.Tinfo;


 
 


