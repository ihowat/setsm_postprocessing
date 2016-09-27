function registerStrips(demdir,gcp,sensor)
% REGISTERSTRIPS register DEM strip data to ground control point cloud
%
% registerStrips(demdir,gcp,sensor) where demdir is the dem directory path
% gcp is a structure with fields x,y,z and optional t. sensor is the sensor
% id string. OUtputs a *_reg.txt file with translation parameters and
% statistics. Uses the interative regression method of Nuth and Kaab, 2011

% Ian Howat, Ohio State University
% Version 1: 27-Sep-2016 09:57:42

fdem=dir([demdir,'/*_dem.tif']);
fdem={fdem.name};
fdem=fdem(:);

if isempty(fdem); return; end

imdate=char(fdem);
imdate=imdate(:,6:13);
imdate=datenum(imdate,'yyyymmdd');

% get integer month of each control point
if isfield(gcp,'t'); [~,gcp.mo] = datevec(gcp.t); end

% loop through DEMs
for i=1:length(fdem)
    
    % make outname
    outfile = [demdir,'/',strrep(fdem{i},'_dem.tif','_reg.txt')];
    
    if exist(outfile,'file'); 
        disp('reg file exists, skipping'); 
          continue
    end
    
    fprintf('processing %d of %d, %s\n',i,length(fdem),[demdir,'/',fdem{i}]);
    
    % get polygon from meta file
    metafile=[demdir,'/',strrep(fdem{i},'_dem.tif','_meta.txt')];
    
    c=textread(metafile,'%s','delimiter','\n'); % read metafile into cellstr
    
    r=find(~cellfun(@isempty,strfind(c,'X:'))); % extract x verts
    xv=str2num(strrep(c{r},'X:',''));
    
    r=find(~cellfun(@isempty,strfind(c,'Y:'))); % extract y verts
    yv=str2num(strrep(c{r},'Y:',''));
    
    clear c r
    
    %% Find overlapping GCPs
    
    % fast rectangle crop gcp data to rectangularized DEM boundarys
    n= gcp.x >= min(xv) & gcp.x <=max(xv) &...
        gcp.y >= min(yv) & gcp.y <=max(yv);
    
    if ~any(n); fprintf('0 overlapping points, skipping\n'); continue; end
    
    % polygon crop
    n=find(n); % reduce index to fast crop selection
    
    nn = inpolygon(gcp.x(n),gcp.y(n),xv,yv); % polygon crop
    
    % skip if too few
    if sum(nn) < 4; 
        fprintf('%d overlapping points, too few,skipping\n',sum(nn)); 
        continue; 
    end
    
    n=n(nn); clear nn % reduce index to fast crop selection
     
    
    if isfield(gcp,'t');
        %% Month selection
        % In order to reduce seasonal variability in height, select gcps within
        % +/- 1 month of the image date. If this results in too few points,
        % expand range incrementally until enough points are found.
        
        [~,mo]=datevec(imdate(i));% get image integer month to select to gcps
        
        % get absolute difference in months between image and gcps accounting
        % for annual wrap-around
        dmo=abs(gcp.mo(n)-mo);
        dmo(dmo > 6) = 12-dmo(dmo > 6);
        
        
        month_range=0;% initial month range +/-
        nn=0; % initialize index
        while sum(nn) < 4; % test if enough gcps
            % expand range by +/- 1 month
            month_range=month_range+1;
            % search  for gcps in ranage
            nn = dmo <= month_range;
        end
        
        n=n(nn); clear nn % crop index to selection
        fprintf('%d points within %d month(s), performing registration\n',...
            length(n),month_range);
    end
    
    %% Registration with subsetting
    % For fitting to control, we can either load the whole image or load
    % subsets of the image in the neighborhood of each control point with
    % a size large enough to allow for shifting within the expected image
    % registration error. The former is fastest for larger numbers of GCPs
    % and or smaller image sizes. We can select by calculating the total
    % number of pixels to be loaded by each method and use the less.
    
    % subset size expect +/- 20m maximum image displacement
    dd=20;
    
    info=imfinfo([demdir,'/',fdem{i}]); % image dimensions info
    if ceil(2.*dd/info.ModelPixelScaleTag(1)+1).^2 * length(n) < info.Width * info.Height %test
        
        % subsetting loads fewer pixels, so we'll use it
        
        % initialize coordinate subset cells
        x=cell(length(n),1);
        y=cell(length(n),1);
        z=cell(length(n),1);
        
        % subsetting loop
        j=1;
        for j=1:length(n);
            
            a= readGeotiff([demdir,'/',fdem{i}],...
                'map_subset',[ gcp.x(n(j)) - dd, gcp.x(n(j)) + dd, ...
                gcp.y(n(j)) - dd, gcp.y(n(j)) + dd ]);
            
            x{j}=a.x;
            y{j}=a.y;
            a.z(a.z <= -100) = NaN;
            z{j} = a.z;
            clear a
        end
        
    else % load whole image
        z= readGeotiff([demdir,'/',fdem{i}]);
        x=z.x;
        y=z.y;
        z=z.z;
        z(z <= -100) = NaN;
    end
    
    % send to registration fx
    [dtrans,dzall] = ...
        registerDEM2LIDAR(x,y,z,gcp.x(n),gcp.y(n),gcp.z(n));
    
    % find points used in registration
    nn = abs(dzall - nanmedian(dzall)) < nanstd(dzall);
    
    clear x y z

    if isempty(dzall); continue; end
    
    ptiles=[50:5:100];
    dzptiles=prctile(abs(dzall),ptiles);

    fid=fopen(outfile,'w');
    fprintf(fid,'DEM Fileame: %s\n',fdem{i});
    j=1; for j=1:length(sensor); fprintf(fid,'Registration Dataset %d Name: %s\n',j,sensor{j});end
    fprintf(fid,'GCP Month Range (+/- months)=%d\n',month_range);
    fprintf(fid,'# GCPs=%d\n',sum(nn));
    fprintf(fid,'Mean Vertical Residual (m)=%.3f\n',nanmean(dzall(nn)));
    fprintf(fid,'Median Vertical Residual (m)=%.3f\n',nanmedian(dzall(nn)));
    fprintf(fid,'Translation Vector (dz,dx,dy)(m)= %.3f, %.3f, %.3f \n',dtrans);
    fprintf(fid,'Vertical Deviation Percentiles(m):\n');
    j=1; for j=1:length(ptiles);
        fprintf(fid,'%dth percentile=',ptiles(j)); fprintf(fid,'%.3f ',dzptiles(j)); fprintf(fid,'\n');
    end
   fclose(fid);
    clear n nn dzall ptiles dzptiles
    fprintf('\n');
end









