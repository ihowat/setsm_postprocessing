function tileRegMeta(f,sensors)
% tileRegMeta: write meta file for ArcticDEM mosaic tile registration stats
%
%   tileRegMeta(f,sensors) where f is either the filename or matfile handle of the
%   tile matfile. Sensors is a cellstr of sensor tags. If sensors ={'neigbor align'}
%   a truncated format will be used.

%% first test if input arg is either a valid filename or mat file handle
if ischar(f) % it's a string, might be a filename
    if exist(f,'file') % yup its a file
        m = matfile(f); % load it
    else
        error('File does not exist');
    end
elseif isvalid(f) % not a string, is it a valid file handle?
    m = f; % yes, so set m to f and get the filename for f
    f = m.Properties.Source;
else
    error('input arg must be a filename or valid matfile handle')
end

% if sensors is passed as a string, change it to cellstr for consistency
if ~iscell(sensors); sensors={sensors}; end

% make outname
outfile=strrep(f,'_dem.mat','.txt');

% test if exists, skip if does
if exist(outfile,'file') 
    fprintf('%s exists, skipping\n',outfile);
    return;
end

%% Alternate Write for Neighbor Align Tile
if length(sensors)==1
    if strcmpi(sensors{1},'Neighbor Align') == 1
        j=1;
        
        %% Create meta file and write universal info
        fid=fopen(outfile,'w');
        fprintf(fid,'DEM Filename: %s\n',f);
        fprintf(fid,'Registration Dataset %d Name: %s\n',j,sensors{j});

        for i=1:length(m.rmse)
            fprintf(fid,'Statistics for Coregistration Cluster %d:\n',i);
            fprintf(fid,'# GCPs=NaN\n');
            fprintf(fid,'Mean Vertical Residual (m)=%.3f\n',m.rmse(:,i));
            fprintf(fid,'Median Vertical Residual (m)=NaN\n');
            fprintf(fid,'Translation Vector (dz,dx,dy)(m)= %.3f, %.3f, %.3f \n',m.dtrans(:,i)');
            fprintf(fid,'Vertical Deviation Percentiles(m):\n');
            fprintf(fid,'NaN\n');
            fprintf(fid,'\n');
        end
    fclose(fid);
    return
    end
end

%% double check to make sure this isnt an unregfile sent here by mistake
if isempty(whos(m,'dzall'))
    error('oops, I think you sent an unregistered tile here by mistake'); 
end

%% Create meta file and write universal info
fid=fopen(outfile,'w');
fprintf(fid,'DEM Filename: %s\n',f);
for j=1:length(sensors)
    fprintf(fid,'Registration Dataset %d Name: %s\n',j,sensors{j});
end    
    
%% Read and format strip info
dzall=m.dzall;
dtrans=m.dtrans;

for i=1:length(dzall)
    fprintf(fid,'Statistics for Coregistration Cluster %d:\n',i);

    % find points used in registration
    nn = abs(dzall{i} - nanmedian(dzall{i})) < nanstd(dzall{i});
    
    % calculate percentiles
    ptiles=50:5:100;
    dzptiles=prctile(abs(dzall{i}),ptiles);

    fprintf(fid,'# GCPs=%d\n',sum(nn));
    fprintf(fid,'Mean Vertical Residual (m)=%.3f\n',nanmean(dzall{i}(nn)));
    fprintf(fid,'Median Vertical Residual (m)=%.3f\n',nanmedian(dzall{i}(nn)));
    fprintf(fid,'Translation Vector (dz,dx,dy)(m)= %.3f, %.3f, %.3f \n',dtrans(:,i)');
    fprintf(fid,'Vertical Deviation Percentiles(m):\n');
    for j=1:length(ptiles)
        fprintf(fid,'%dth percentile=',ptiles(j)); fprintf(fid,'%.3f ',dzptiles(j)); fprintf(fid,'\n');
    end
   
    fprintf(fid,'\n');
        
end
 fclose(fid);

