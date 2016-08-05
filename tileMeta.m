function tileMeta(f)
% tileMeta: write meta file for ArcticDEM mosaic tile
%
%   tileMeta(f) where f is either the filename or matfile handle of the
%   tile matfile.

%% first test if input arg is either a valid filename or mat file handle
if isstr(f) % it's a string, might be a filename
    if exist(f,'file') % yup its a file
        m = matfile(f); % load it
    else
        error('File does not exist');
    end
elseif isvalid(f) % not a string, is it a valid file handle?
    m = f; % yes, so set m to f and get the filename for f
    f = m.Properties.Source;
else error('input arg must be a filename or valid matfile handle')
end

% make outname
outfile=strrep(f,'.mat','_meta.txt');

% test if exists, skip if does
if exist(outfile,'file'); 
    fprintf('%s exists, skipping\n',f);
    return;
end
    
%% Read and format strip info
strip.name=m.f;
strip.name=strip.name(:);

% strip out path from file names
for i=1:length(strip.name);
    [~,strip.name{i}]=fileparts(strip.name{i});
end
% remove meta suffix from filenames
strip.name=strrep(strip.name,'_meta','');

strip.dtrans=m.dtrans;
strip.dtrans=strip.dtrans';

strip.rmse=m.rmse;
strip.rmse=strip.rmse(:);

% find nan dtrans, meaning strip wasn't used, and remove from lists
n = any(~isnan(strip.dtrans'))';

strip = structfun(@(x) ( x(n,:) ), strip, 'UniformOutput', false);

%% Get tile version
tileVersion='1.0';
if ~isempty(whos(m,'version')); tileVersion=m.version; end;

%% Write File
fid=fopen(outfile,'w');
fprintf(fid,'ArcticDEM Mosaic Tile Metadata \n');
fprintf(fid,'Creation Date: %s\n',datestr(now));
fprintf(fid,'Version: %s\n',tileVersion);
fprintf(fid,'\n');
fprintf(fid,'Mosaicking Alignment Statistics (meters) in rank order\n');
fprintf(fid,'strip, rmse, dz, dx, dy\n');
i=1;
for i=1:length(strip.name)
    fprintf(fid,'%s %.2f %.4f %.4f %.4f\n',...
        strip.name{i},strip.rmse(i),strip.dtrans(i,:)); 
end
fclose(fid);
    
