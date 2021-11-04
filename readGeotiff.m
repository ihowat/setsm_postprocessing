function I= readGeotiff(name,varargin)
% readGeotiff: read geotiff using imread and assign map info from infinfo.
%
% I = readGeotiff('filename'); reads entire .tif into map structure I
% I = GEOTIFF_READ(...,'pixel_subset', [r0 r1 c0 c1]) reads the subset of
% the image defined by row and column range.
% I = GEOTIFF_READ(...,'map_subset', [x0 x1 y0 y1]); reads the subset of
% the image defined by the map coordinate range.
% output:
% I.z, image data
% I.x, x coordinate in map
% I.y, y coordinate in map
% I.info, misc. info
% imshow(I.z, 'xdata', I.x, 'ydata', I.y);
% shows image with map coordinate
% Version by Yushin Ahn, ahn.74@osu.edu
% Glacier Dynamics Laboratory, 
% Byrd Polar Resear Center, Ohio State University 
% Referenced enviread.m (Ian Howat)

n = find(strcmpi(varargin,'reportMapSubsetError'));
if ~isempty(n)
    reportMapSubsetError = true;
else
    reportMapSubsetError = false;
end

Tinfo       = imfinfo(name);
info.cols   = Tinfo.Width;
info.rows   = Tinfo.Height;
info.imsize = Tinfo.Offset;
info.bands  = Tinfo.SamplesPerPixel;

% pull out map fields in case overviews make them structures. 
ModelPixelScaleTag = Tinfo.ModelPixelScaleTag;
ModelTiepointTag = Tinfo.ModelTiepointTag;
info.map_info.dx    = ModelPixelScaleTag(1);
info.map_info.dy    = ModelPixelScaleTag(2);
info.map_info.mapx  = ModelTiepointTag(4);
info.map_info.mapy  = ModelTiepointTag(5);

subrows = [1 info.rows];
subcols = [1 info.cols];

minx = info.map_info.mapx;
maxy = info.map_info.mapy;

x = minx + ((0:info.cols-1).*info.map_info.dx);
y = maxy - ((0:info.rows  -1).*info.map_info.dy);

mapInfoOnlyFlag=0;

% get subset bounds if specified
for i=1:length(varargin);
	    
    if strcmpi(varargin{i},'pixel_subset');
        i=i+1;
	subrows = varargin{i}(3:4);
        subcols = varargin{i}(1:2);
        
    elseif strcmpi(varargin{i},'map_subset');
        i=i+1;
        map_subset  = varargin{i};
        subcols = (map_subset(1:2)-info.map_info.mapx)./...
            info.map_info.dx+1;
        subrows = (info.map_info.mapy - map_subset([4,3]))./...
            info.map_info.dy+1;
        subcols= round(subcols);
        subrows = round(subrows);

        if reportMapSubsetError
            if any(subcols < 1) || any(subrows < 1) || any(subcols > info.cols) || any(subrows > info.rows)
                map_subset
                error('readGeotiff map_subset extends beyond raster extent (x=[%d,%d], y=[%d,%d]): %s',...
                    x(1),x(end),y(end),y(1),name)
            end
        end
        subcols(subcols < 1) = 1;
        subrows(subrows < 1) = 1;
        subcols(subcols > info.cols) = info.cols;
        subrows(subrows > info.rows) = info.rows;
        
     elseif strcmpi(varargin{i},'mapinfoonly');
         mapInfoOnlyFlag=1;
    end

end

% define map vectors
I.x = x(subcols(1):subcols(2));
I.y = y(subrows(1):subrows(2));
I.z = [];

if ~mapInfoOnlyFlag
    % read image
    I.z = imread(name,'PixelRegion',{subrows,subcols});
%    try
%        I.z = imread(name,'PixelRegion',{subrows,subcols});
%    catch ME
%        [stripdir, stripfname, stripfext] = fileparts(name);
%        [~, stripdname, ~] = fileparts(stripdir);
%        altname = fullfile('/mnt/pgc/data/scratch/erik/s2s_reruns/strips/2m', stripdname, [stripfname, stripfext]);
%        if isfile(altname)
%            fprintf('Switching to alternate file location: %s\n', altname)
%            I.z = imread(altname,'PixelRegion',{subrows,subcols});
%        else
%            rethrow(ME)
%        end
%    end
end


I.info = info;
I.Tinfo = Tinfo;





