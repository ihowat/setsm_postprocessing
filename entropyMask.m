function M = entropyMask(orthoFile)
% entropyMask classify areas of low entropy in an image such as water
%
% M = entropyMask(orthoFile) returns the low entropy classification mask
% from the geotif image in orthoFile. Also checks whether wvc was applied
% from the metafile.
%
% Ian Howat,ihowat@gmail.com, Ohio State
% 13-Apr-2017 10:41:41

% parameters
pres = 8; % resample to this resolution for processing for speed and smooth
Jthresh = 0.2; % minimum entropy threshold. 0.2 seems to be good for water
minPix = 1000; % clusters of mask and void pixels in the resampled
               % image less than this will be removed

% Get meta file info
metaFile = strrep(orthoFile,'ortho.tif','meta.txt');

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

% read ortho
or=readGeotiff(orthoFile);
res = or.x(2) - or.x(1);
or = or.z;
sz = size(or);

bg = or == 0; % image background mask

% resize ortho to pres
if res ~= pres
    or = imresize(or,res./pres);
end

% subtraction image
or_subtraction =  movmax(or,5) - movmin(or,5);
 
if ~wvcFlag; or_subtraction = uint8(or_subtraction); end

 % entropy image
J = entropyfilt(or_subtraction,ones(5));

M = J < Jthresh;

M = bwareaopen(M,minPix);
M = ~bwareaopen(~M,minPix);

% resize ortho to 8m
if res ~= pres
    M = imresize(M,sz,'nearest');
end

M(bg) = false;

