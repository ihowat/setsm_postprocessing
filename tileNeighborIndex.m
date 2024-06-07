function tileNeighborIndexFile=tileNeighborIndex(tileDir,varargin)

if ~isfolder(tileDir)
    error("input 'tileDir' folder does not exist: %s", tileDir);
end

org = 'osu';
resolution = '10m';
outfile = [];
priority_suffix = '_reg.mat';
secondary_suffix = '.mat';

if ~isempty(varargin)
    org_choices = {'osu', 'pgc'};
    n = find(strcmpi('org', varargin));
    if ~isempty(n)
        org = varargin{n+1};
    end
    if ~ismember(org, org_choices)
        error("'org' must be one of the following, but was '%s': {'%s'}", org, strjoin(org_choices, "', '"));
    end

    resolution_choices = {'10m', '2m'};
    n = find(strcmpi('resolution', varargin));
    if ~isempty(n)
        resolution = varargin{n+1};
    end
    if ~ismember(resolution, resolution_choices)
        error("'resolution' must be one of the following, but was '%s': {'%s'}", resolution, strjoin(resolution_choices, "', '"));
    end

    n = find(strcmpi('outfile', varargin));
    if ~isempty(n)
        outfile = varargin{n+1};
    end

    n = find(strcmpi('priority_suffix', varargin));
    if ~isempty(n)
        priority_suffix = varargin{n+1};
    end

    n = find(strcmpi('secondary_suffix', varargin));
    if ~isempty(n)
        secondary_suffix = varargin{n+1};
    end
end

if strcmp(org, 'osu')
    use_res_in_outname = false;
else
    use_res_in_outname = true;
end
if strcmp(resolution, '2m')
    use_res_in_outname = true;
end

fileNames = selectTilesForMerging(tileDir,varargin{:});

% find neighbors
[ntop,nright,ntopRight,nbottomRight] = findNeighborTiles(fileNames);

nN = nan(length(fileNames),8);
nN(ntop(:,1),1) = ntop(:,2); %top
nN(ntop(:,2),2) = ntop(:,1); % bottom
nN(nright(:,2),3) = nright(:,1); %left
nN(nright(:,1),4) = nright(:,2); %right
nN(nbottomRight(:,2),5) = nbottomRight(:,1); %top-left
nN(ntopRight(:,1),6) = ntopRight(:,2); %top-right
nN(ntopRight(:,2),7) = ntopRight(:,1); %bottom-left
nN(nbottomRight(:,1),8) = nbottomRight(:,2); %bottom-right

if ~isempty(outfile)
    tileNeighborIndexFile=outfile;
elseif use_res_in_outname
    tileNeighborIndexFile=[tileDir,'/tileNeighborIndex_',resolution,'.mat'];
else
    tileNeighborIndexFile=[tileDir,'/tileNeighborIndex.mat'];
end
fprintf('writing %s\n', tileNeighborIndexFile)
save(tileNeighborIndexFile,'fileNames','nN')
