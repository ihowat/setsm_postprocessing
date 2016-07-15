function enviwrite(varargin)
% ENVIWRITE write an array to an ENVI binary raster file and header.
%   ENVIWRITE('filename',Z) writes the array Z to the ENVI default binary
%   format. Also written is the ascii ENVI header file 'filename.hdr'.
%   ENVIWRITE('filename',X,Y,Z,'proj','projection name') uses 

%enviwrite(file,x,y,z,fmt);
%creates envi ascii data and header files from cartesian x,y,z arrays.
%keyboard
minarg = 4;
maxarg = 12;

% Check input
error(nargchk(minarg,maxarg,nargin));

% requisite arguments
file = varargin{1};
x = varargin{2};
y = varargin{3};
z = varargin{4};
for i = 1:size(z,3);
    zr(:,:,i) = z(:,:,i)';
end

z = zr;

% Default option values
fmt = 4;                 % floating point format
proj = 'UTM';
hemi = 'north';
datum = 'wgs-84';
units = 'Meters';
proj_info = [];
% Alteration of default values
for n = (minarg+1:2:nargin-1);
    % Checks that every value is a numeric value and every property
    % is a sting
    if ~isstr(varargin{n})
        error('Invalid parameter/value pair arguments.');
        % These properties take string arguments otherwise it must be
        % numeric value
    elseif ~isnumeric(varargin{n+1})    &&...
           ~strcmpi(varargin{n},'proj') &&...
           ~strcmpi(varargin{n},'hemi') &&...
           ~strcmpi(varargin{n},'datum')&&...
           ~strcmpi(varargin{n},'units')
   
     msg = sprintf(['Bad value for property: '...
            '''' varargin{n} '''\n unknown option.']);
        error(msg);
    end

    % Sets property value
    switch lower(varargin{n})
        case 'format'
            fmt = varargin{n+1};
        case 'zone'
            zone = varargin{n+1};
        case 'proj'
            proj = varargin{n+1};
        case 'hemi'
            hemi = varargin{n+1};
        case 'datum'
            datum = varargin{n+1};
        case 'units'
            units = varargin{n+1};
    end
end


switch lower(proj)
    
    case 'polar stereo north';
        proj_info = '31, 6378137.0, 6356752.3, 70.000000, -45.000000, 0, 0, WGS-84, Polar Stereo North, units=Meters';

     case 'polar stereo south';   
        proj_info = '31, 637813.0, 635674.5, -70.000000, 0.000000, 0, 0, WGS-84, Polar Stereo South, units=Meters';
        
    case 'bamber'
         proj_info = '31, 6378137.0, 6356752.3, 70.000000, -39.000000, 1, 1, WGS-84, bamber, units=Meters';
    case 'bamber2'
         proj_info = '31, 6378137.0, 6356752.3, 71.000000, -39.000000, 1, 1, WGS-84, bamber, units=Meters';
end

if size(x,1) > 1 && size(x,2) > 1;
    x = x(1,:);
else
    x = x(:);
end;

if size(y,1) > 1 && size(y,2) > 1;
    y = y(:,1);
else
    y = y(:);
end;


d = [datestr(now,8),' ',datestr(now,3),' ',datestr(now,7),' ',datestr(now,13),' ',datestr(now,10)];
[samples,lines,bands] = size(z);
x1 = min(x);
y1 = max(y);
dx = abs(diff(x(1:2)));
dy = abs(diff(y(1:2)));


switch fmt
    case 1; fmtstr = 'uint8';
    case 2; fmtstr = 'int16';
    case 3; fmtstr = 'int32';
    case 4; fmtstr = 'single';
    case 5; fmtstr = 'double';
    case 6; error('cant write complex numbers') %complex single
    case 9; error('cant write complex numbers') %complex double
    case 11; fmtstr = 'int8'; fmt = 1; %envi can't do signed byte
    case 12; fmtstr = 'uint16';
    case 13; fmtstr = 'uint32';
    case 14; fmtstr = 'int64';
    case 15; fmtstr = 'uint64';
end
         
fid = fopen(file,'w');
count = fwrite(fid,z,fmtstr); fclose all;


fid = fopen([file,'.hdr'],'w');%
fprintf(fid,'%s\n','ENVI'); %
fprintf(fid,'%s\n','description = {');
fprintf(fid,'%s\n',['  Exported Matlab Array [',d,']}']);
fprintf(fid,'%s\n',['samples = ',num2str(samples)]);
fprintf(fid,'%s\n',['lines   = ',num2str(lines)]);
fprintf(fid,'%s\n',['bands   = ',num2str(size(z,3))]);
fprintf(fid,'%s\n','header offset = 0');
fprintf(fid,'%s\n','file type = ENVI Standard');
fprintf(fid,'%s\n',['data type = ',num2str(fmt)]);
fprintf(fid,'%s\n','interleave = bsq');
fprintf(fid,'%s\n','sensor type = Unknown');
fprintf(fid,'%s\n','byte order = 0');
if strcmpi(proj,'UTM'); 
fprintf(fid,'%s\n',['map info = { ',proj,', 1.000, 1.000, ',num2str(x1),', ',num2str(y1),', ',num2str(dx),', ',num2str(dy),', ',num2str(zone),', ',hemi,', ',datum,', units=',units,'}']);
else
    fprintf(fid,'%s\n',['map info = { ',proj,', 1.000, 1.000, ',num2str(x1),', ',num2str(y1),', ',num2str(dx),', ',num2str(dy),',',datum,', units=',units,'}']);
    fprintf(fid,'%s\n',['projection info = {',proj_info,'}']);
end;
fprintf(fid,'%s\n','wavelength units = Unknown');
fclose(fid);% close the file
