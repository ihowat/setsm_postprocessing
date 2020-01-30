function tileMetav4(f)
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

fileAtts= dir(m.Properties.Source);

% test if exists, skip if does
if exist(outfile,'file') 
    fprintf('%s exists, skipping\n',outfile);
    return;
end
    
%% Read and format strip info
strip.name=unique(m.stripList);
strip.name=strip.name(:);

%% Get tile version
tileVersion='4.0';
if ~isempty(whos(m,'version')); tileVersion=m.version; end

fprintf('writing meta\n')

%% Write File
fid=fopen(outfile,'w');
fprintf(fid,'ArcticDEM Mosaic Tile Metadata \n');
fprintf(fid,'Creation Date: %s\n',fileAtts.date);
fprintf(fid,'Version: %s\n',tileVersion);
fprintf(fid,'\n');
fprintf(fid,'ICESat N points: %d\n',m.icesat_dz_N);
fprintf(fid,'ICESat median offset: %.2f\n',m.icesat_dz_med);
fprintf(fid,'ICESat median absolute deviation: %.2f\n',m.icesat_dz_mad);
fprintf(fid,'\n');
fprintf(fid,'List of strips used in mosaic:\n');
i=1;
for i=1:length(strip.name)
    fprintf(fid,'%s\n',...
        strip.name{i}); 
end
fclose(fid);
    
