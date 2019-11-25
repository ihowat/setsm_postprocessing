function info = read_envihdr(hdrfile)
% READ_ENVIHDR read and return ENVI image file header information.
%   INFO = READ_ENVIHDR('HDR_FILE') reads the ASCII ENVI-generated image
%   header file and returns all the information in a structure of
%   parameters. 
%
% Output: 
%  * Info - struct with fields provided in the ENVI file. ENVI header 
%             format requires the following fields:
%    * samples - number of samples in the image (columns)
%    * lines - number of lines in the image (rows)
%    * bands - number of bands in the image. If all 3 dimensions are
%      provided than info.size will be created holding 
%      [info.lines info.samples info.bands]
%    * data_type - data type of the image stored as an integer in 1-15 
%      range. If provided than info.format will be created holding string 
%      with Matlab's type name.  
%    * interleave -  file band interleave type; either bip, bsq, or bil are 
%      possible
%    * byte_order - byte order (0 is little endian [least significant
%       byte first], 1 is big endian [most significant byte first]). 
%      If provided than info.machine will be created holding either 
%      'ieee-le' or 'ieee-be' string.
%
%   Example 1:
%   >> info = read_envihdr('my_envi_image.hdr')
%   info =
%          description: [1x101 char]
%              samples: 658
%                lines: 749
%                bands: 3
%                 size: [749 658 3]
%        header_offset: 0
%            file_type: 'ENVI Standard'
%            data_type: 4
%              format : 'single'
%           interleave: 'bsq'
%          sensor_type: 'Unknown'
%           byte_order: 0
%             map_info: [1x1 struct]
%      projection_info: [1x102 char]
%     wavelength_units: 'Unknown'
%           pixel_size: [1x1 struct]
%           band_names: [1x154 char]
% Example 2:
% >> info = read_envihdr('my_envi_image.hsi');
% >> Z = multibandread(gFile, info.size, [info.format '=>double'], ...
% >>           info.header_offset, info.interleave, info.machine);
%
% Author: Jarek Tuszynski (jaroslaw.w.tuszynski@saic.com)
% License: BSD (Berkeley Software Distribution)
%
% See Also
% * MATLAB function multibandread
% * http://geol.hu/data/online_help/ENVI_Header_Format.html

%% If file does not have 'hdr' extension than check if there is an matching 
% header file
[FilePath, FileRoot, FileExt] = fileparts(hdrfile);
if ~strcmp(FileExt, '.hdr')
  fname = hdrfile;
  hdrfile = [fname '.hdr']; % add '.hdr' to the file name
  if ~exist(hdrfile,'file')
    hdrfile = [FilePath '\' FileRoot '.hdr']; % replace extension
    if ~exist(hdrfile,'file')
      hdrfile = fname;
    end
  end
end

%% Load whole header file into a string
fid = fopen(hdrfile);
if fid<0, error('%s does not exist. \n', hdrfile); end
str = fread(fid,'uint8=>char')';
fclose(fid);

%% split string into lines 
flag = 0;
str(str==10) = 13;
str=strrep(str,char([13 13]), char(13));
for i = 1:length(str)
  switch str(i)
    case '{'
      flag=flag+1;
    case 13
      if (flag), str(i)=10; end 
    case '}'
      flag=flag-1;
  end
end
lines = textscan(str,'%s','Delimiter',char(13));
lines = lines{1};

%% parse each line into a field of a struc
info = [];
for iLine=1:length(lines)
  [info param] = ParseLine(lines{iLine}, info);
  if ~isempty(param) && ischar(info.(param)) && nnz(info.(param)=='{')
    % if "{" found than parse one more level
    line = info.(param);
    if nnz(info.(param)=='=')==0 % string has no "=" -> check if it is numeric array
      line(line<32) = [];
      line(line == '{') = '[';
      line(line == '}') = ']';
      num = str2num(line); %#ok<ST2NM>
      if isnumeric(num) && ~isempty(num)
        info.(param) = num;
      elseif nnz(info.(param)==',')>0 % string has "," -> split into string cell array
        line(line == '[' | line == ']') = [];
        lines2 = textscan(line,'%s','Delimiter',',');
        info.(param) = lines2{1};
      end
    else                         % string has "="
      line(line == '{' | line == '}') = [];
      lines2 = textscan(line,'%s','Delimiter',char(10));
      lines2 = lines2{1};
      info2 = [];
      for jLine=1:length(lines2)
        info2 = ParseLine(lines2{jLine}, info2);
      end
      info.(param) = info2;
    end
  end
end

%% create info.size
if isfield(info, 'lines') && isfield(info, 'samples') && isfield(info, 'bands')
  info.size = [info.lines info.samples info.bands];
end

%% fix 'byte_order' field
if isfield(info,'byte_order')
  switch info.byte_order
    case 0
      info.machine = 'ieee-le';
    case 1
      info.machine = 'ieee-be';
    otherwise
      info.machine = 'n';
  end
end

%% fix 'data_type' field
if isfield(info,'data_type')
  info.iscomplex=false; %if it is complex
  switch info.data_type
    case 1
      info.format = 'uint8';
    case 2
      info.format= 'int16';
    case 3
      info.format= 'int32';
    case 4
      info.format= 'single';
    case 5
      info.format= 'double';
    case 6
      info.iscomplex=true;
      info.format= 'single';
    case 9
      info.iscomplex=true;
      info.format= 'double';
    case 12
      info.format= 'uint16';
    case 13
      info.format= 'uint32';
    case 14
      info.format= 'int64';
    case 15
      info.format= 'uint64';
    otherwise
      error(['File type number: ',num2str(dtype),' not supported']);
  end
end


function [struc param] = ParseLine(line, struc)
  param='';
  eqsn = find(line=='=',1,'first'); % find = sign
  if ~isempty(eqsn)
    param = strtrim(line(1:eqsn-1));
    param(strfind(param,' ')) = '_';
    param = genvarname(param);         % split line into param
    value = strtrim(line(eqsn+1:end)); % and value
    if strcmp(param,upper(param))      % is all letters are upper case than ... 
      param = lower(param);            % convert to lower case string
    end

    % save values as fields of a struct
    try
      struc.(param) = eval(value);    % try converting to numbers, etc.
    catch ME %#ok<NASGU>
      struc.(param) = value;
    end
  else
    struc.CONTENT = line;
  end