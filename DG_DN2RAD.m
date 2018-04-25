function [L, effectiveBandwith, abscalFactor, gain, offset] = DG_DN2RAD(DN,varargin) 
% DG_DN2RAD converts DG DN images to top-of-atmosphere radiance
%
% L =  = DG_DN2RAD(DN, satID, effectiveBandwith, abscalFactor) applies the
% conversion using the supplied factors with a table look-up for the satID.
% The output L is top-of-atmosphere radiance in units of Wµm^-1 m^-2 sr^-1.
%
% L =  = DG_DN2RAD(DN,xmlFile) reads the factors from the supplied xml file
%
% [L, effectiveBandwith, abscalFactor, gain, offset] = DG_DN2RAD(...)
% returns the scaling parameters used.
 

% parse inputs
switch length(varargin)
    case 1
        xmlFile=varargin{1};
        [satID,effectiveBandwith,abscalFactor] = getParams(xmlFile);
        if strcmp(satID,'QB2'); satID = 'QB02'; end
        if strcmp(satID,'IKO'); satID = 'IK01'; end
    case 3
        satID              = varargin{1};
        effectiveBandwith  = varargin{2};
        abscalFactor       = varargin{3};
    otherwise
        error('must be 2 or 4 input arguments')
end
     
% Values from:
% https://dg-cms-uploads-production.s3.amazonaws.com/uploads/document/file/209/DGConstellationAbsRadCalAdjustmentFactors_2015v2.pdf
sensor = {'WV03', 'WV02', 'GE01', 'QB02', 'IK01', 'WV01'};
gain = [  0.923	0.96	0.978	0.876	0.907	1.016];
offset = [ -1.7	-2.957	-1.948	-2.157	-4.461	-3.932];

n = find(strcmp(satID,sensor));
gain = gain(n);
offset = offset(n);

%convert to single
DN = single(DN);

% set nodata to nan
DN(DN == 0) = NaN;

% calculate radiance
L = gain.*DN.*(abscalFactor./effectiveBandwith) + offset;


function [satID,effectiveBandwith, abscalFactor] = getParams(xmlFile)

paramstr='EFFECTIVEBANDWIDTH';
val = readFromXml(xmlFile,paramstr);
effectiveBandwith = str2double(val);

paramstr='ABSCALFACTOR';
val = readFromXml(xmlFile,paramstr);
abscalFactor  = str2double(val);

paramstr='SATID';
satID = readFromXml(xmlFile,paramstr);

function val = readFromXml(xmlFile,paramstr)

c=textread(xmlFile,'%s','delimiter','\n');
r=find(~cellfun(@isempty,strfind(c,['<',paramstr,'>'])));
val=deblank(strrep(c{r(1)},['<',paramstr,'>'],''));
val=deblank(strrep(val,['</',paramstr,'>'],''));

