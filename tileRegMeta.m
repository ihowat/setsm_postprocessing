function tileRegMeta(f)
% tileRegMeta: write meta file for ArcticDEM mosaic tile registration stats
%
%   tileRegMeta(f) where f is either the filename or matfile handle of the
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
outfile=strrep(f,'.mat','.txt');

% test if exists, skip if does
if exist(outfile,'file'); 
    fprintf('%s exists, skipping\n',outfile);
    return;
end

%% double check to make sure this isnt an unregfile sent here by mistake
if isempty(whos(m,'dzall')); error('oops, I think you sent an unregistered tile here by mistake'); end

sensors={'GLA14_rel34'};


%% Create meta file and write universal info
fid=fopen(outfile,'w');
fprintf(fid,'DEM Filename: %s\n',f);
for j=1:length(sensors); fprintf(fid,'Registration Dataset %d Name: %s\n',j,sensors{j});end    
    
%% Read and format strip info
dzall=m.dzall;
dtrans=m.dtrans;

for i=1:length(dzall)
    fprintf(fid,'Statistics for Coregistration Cluster %d:\n',i);

    % find points used in registration
    nn = abs(dzall{i} - nanmedian(dzall{i})) < nanstd(dzall{i});
    
    % calculate percentiles
    ptiles=[50:5:100];
    dzptiles=prctile(abs(dzall{i}),ptiles);

    fprintf(fid,'# GCPs=%d\n',sum(nn));
    fprintf(fid,'Mean Vertical Residual (m)=%.3f\n',nanmean(dzall{i}(nn)));
    fprintf(fid,'Median Vertical Residual (m)=%.3f\n',nanmedian(dzall{i}(nn)));
    fprintf(fid,'Translation Vector (dz,dx,dy)(m)= %.3f, %.3f, %.3f \n',dtrans(:,i)');
    fprintf(fid,'Vertical Deviation Percentiles(m):\n');
    for j=1:length(ptiles);
        fprintf(fid,'%dth percentile=',ptiles(j)); fprintf(fid,'%.3f ',dzptiles(j)); fprintf(fid,'\n');
    end
   
    fprintf(fid,'\n');
        
end
 fclose(fid);

