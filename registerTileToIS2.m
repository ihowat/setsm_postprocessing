function registerTileToIS2(tileFile,is2TileFile)


m=matfile(tileFile);

% standard resizeFraction is for 10m, scale for 2m and other res
resizeFraction_10m = 1.0;
res = m.x(1,2) - m.x(1,1);
resizeFraction = min(1.0, resizeFraction_10m * (res/10));

% Erik decided to disable the automatic downsample capability for the
% 2022 REMAv2 mosaic production after a brief analysis showed that
% downsampling the 2m during registration resulted in an inconsitent shift
% between the 10m and 2m versions of mosaic tiles.
resizeFraction = 1.0;


% make sure data array altering post-process steps have not been applied
m_struct = m;
m_varlist = who(m_struct);

dzfitApplied = ismember('dzfitApplied', m_varlist) && m_struct.dzfitApplied;
adjusted = ismember('adjusted', m_varlist) && m_struct.adjusted;

mergedTop = ismember('mergedTop', m_varlist) && m_struct.mergedTop;
mergedBottom = ismember('mergedBottom', m_varlist) && m_struct.mergedBottom;
mergedLeft = ismember('mergedLeft', m_varlist) && m_struct.mergedLeft;
mergedRight = ismember('mergedRight', m_varlist) && m_struct.mergedRight;

if dzfitApplied
    error('unregistered tile somehow has dzfit applied, please investigate')
end
if adjusted
    error('unregistered tile has adjust flag true, cannot register')
end
if mergedTop || mergedBottom || mergedLeft || mergedRight
    error('unregistered tile has had at least one side merged, cannot register')
end


if any(strcmp(fields(m),'unreg'))
    fprintf('unreg field already exists, clearing: %s\n', tileFile)
    m.Properties.Writable = true;
    m.unreg = [];
end
if any(strcmp(fields(m),'reg'))
    fprintf('reg field already exists, clearing: %s\n', tileFile)
    m.Properties.Writable = true;
    m.reg = [];
end

% skip if is2 file doesnt exist for this tile
if ~exist(is2TileFile,'file')
    fprintf('%s doesnt exist, skipping\n',is2TileFile)
  return
end

% Read altimetry file
is2=load(is2TileFile,'x','y','z','t');

try
    %load DEM
    load(tileFile,'x','y','z','land');
catch
    warning('cant read %s\n',tileFile)
    return
end

% remove is2 points that fall outside of the (2m) tile
n = is2.x >= min(x) & is2.x <= max(x) & is2.y >= min(y) & is2.y <= max(y);
is2 = structfun( @(x) x(n), is2, 'uniformoutput', 0);

if resizeFraction ~= 1.0
    x = imresize(x, resizeFraction);
    y = imresize(y, resizeFraction);
    z = imresize(z, resizeFraction);
    land = imresize(land, resizeFraction, 'nearest');
end

if any(~land(:))
    landFlag = interp2(x,y,land,is2.x,is2.y,'*nearest',0);
    landFlag=landFlag==1;
    is2 = structfun( @(x) x(landFlag), is2, 'uniformoutput', 0);
end

% get dem elevations at is2 points
is2.dem = interp2(x,y,z,is2.x,is2.y,'*linear');

% pre-registration residuals
dz0 = is2.dem - is2.z;

% pre-registrations stats
unreg.N = sum(~isnan(dz0));
unreg.mn = nanmean(dz0);
unreg.md = nanmedian(dz0);
unreg.sd = nanstd(dz0);
unreg.le68 = prctile(abs(dz0),68);
unreg.le90 = prctile(abs(dz0),90);

m=matfile(tileFile);
m.Properties.Writable = true;
m.unreg=unreg;

% registration
[p,dzr,~,sigma_dtrans,p68,p90] = registerDEM2GCP(x,y,z,is2.x,is2.y,is2.z);

% check if reg failed
if ~isempty(sigma_dtrans)
    % check of horizontal registation failed and subtract median if yes
    if p(1) == 0
        p(1) = nanmedian(dzr);
        % set error to standard erro of mean
        perr=[nanstd(dzr)./sqrt(length(dzr(~isnan(dzr))));0;0];
        
        % subtract median from residual and get LE's
        dzr = dzr - nanmedian(dzr);
        p68=prctile(abs(dzr(~isnan(dzr))),68);
        p90=prctile(abs(dzr(~isnan(dzr))),90);
    else
        perr = sqrt(sum(sigma_dtrans.^2,2));
    end
    
    reg.p = p(:)';
    reg.perr = perr(:)';
    reg.N = sum(~isnan(dzr));
    reg.mn = nanmean(dzr);
    reg.md = nanmedian(dzr);
    reg.sd = nanstd(dzr);
    reg.le68 = p68(end);
    reg.le90 = p90(end);
    reg.tavg = nanmean(is2.t(~isnan(dzr)));
    
    m.reg = reg;

    m.aux_regIS2_is2TileFile = is2TileFile;
end




