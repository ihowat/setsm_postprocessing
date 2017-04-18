function strips2tile(varargin)
% STRIPS2TILE mosaic strips to rectangular tile
%
%[x,y,z,mt,or,dy,f,dtrans,rmse] = ...
%           strips2tile(meta,tilex0, tilex1, tiley0, tiley1,res,outname)
%
% where meta is the database structure, tilex and tiley are the coordinate
% ranges of the tile, buff is an edge buffer to add  to these ranges in pixels
% res is the tile resolution in m and outname is the tile output file name.
%
%[...] = strips2tile(...,'options') where the following options are
%availble:
%       'disableReg' disables use of any existing GCP registration
%       files/data
%       'disableCoregTest' disables the coregistraiton optimization test
%       to save time
%       'mergeMethod','methodstring', specifies the mergeing method as
%       'feather' (default), 'underprint' or 'warp'. See addStrip2Mosaic
%       for descriptions.
%       'mergeMethodReg','methodstring', specifies the mergeing method for
%       registered strips.'feather' (default), 'underprint' or 'warp'.
%       'buffer' number of edge buffer pixels to add to the range, default=100
%
% subfunctions: stripSearch,regStrips2Tile,initializeTile,readreg,addStrip2Mosaic
%
% Procedure description: Strips are added to the mosaic to form
% 'coregistration clusters'. Strips with a priori registration are assigned
% a cluster number of 1, as are any strips that are coregistered to that
% cluster. Strips without a priori registration and that cannot be
% coregistered to cluster 1 will be added as a new cluster, starting at #2.
% Thus each cluster # represents a group of coregistered strips, with C=1
% being absolute (reference to a priori control), and C>=2 being relative
% (internally coregistered). Therefore, any subsequent transformation
% shoudl be applied individually to each cluster.
%
%   Ian Howat, Ohio State University
%   Version 3.0; 10-Feb-2017 13:56:01

tileVersion='3.0';

%% Set Parameters & Parse args
% set Y2K as day 0 for the daynumber grid
dy0 = datenum('jan 1 2000');

% option to disable using control registration files
disableReg=false;

% option to disable coregistration optimization test
disableCoregTest=false;

% default merging method
mergeMethod='feather';

%default edge buffer
buff=100;

% default merging method for registered files
mergeMethodReg='underprint';

if nargin >= 7 % number of argins needed for creating a new mosaic from scratch
    
    meta    = varargin{1};
    tilex0  = varargin{2};
    tilex1  = varargin{3};
    tiley0  = varargin{4};
    tiley1  = varargin{5};
    res     = varargin{6};
    outname = varargin{7};
    
    % test varargin for flags
    if ~isempty(varargin)
        if any(strcmpi('disableReg',varargin))
            fprintf('GCP registration disabled\n')
            disableReg=true;
        end
        
        if any(strcmpi('disableCoregTest',varargin))
            fprintf('Coregistration testing disabled\n')
            disableCoregTest=true;
        end

        if any(strcmpi('mergeMethod',varargin))
            n=find(strcmpi('mergeMethod',varargin));
            mergeMethod=varargin{n+1};
        end
        
        if any(strcmpi('mergeMethodReg',varargin))
            n=find(strcmpi('mergeMethodReg',varargin));
            mergeMethodReg=varargin{n+1};
        end
        
        if any(strcmpi('buffer',varargin))
            n=find(strcmpi('buffer',varargin));
            buff=varargin{n+1};
        end
        
    end
    
else
    error('incorrect number of input arguments')
end

fprintf('Using registered strip merge method: %s\n',mergeMethodReg)
fprintf('Using floating strip merge method: %s\n',mergeMethod)

% Filter strips not overlapping this tile
meta = stripSearch(meta,tilex0,tilex1,tiley0,tiley1);

% if all meta cleared, return
if isempty(meta); return; end

% Initialize/Restart Mosaic
[x,y,c,C,N,m,meta] = initializeTile(...
    meta,tilex0,tilex1,tiley0,tiley1,buff,res,outname,tileVersion,disableReg,...
    disableCoregTest,mergeMethod,mergeMethodReg);

% Add Registered Strips
if ~disableReg % check for disable registration argument
    [meta,m,N,C] = regStrips2Tile(meta,m,N,C,mergeMethodReg,dy0);
end

% Build Grid Point Index Field
meta = buildGridPointInd(meta,x,y,N);

%% Sequential Floating Strip Addition Loop
% Sequentially add strips to the mosaic by selecting the file with the most
% new data coverage, with enough overlap to coregister to exisiting data.

minNewPixels = 1000/res;
minOverlapPixels = 1000/res;

% add dummy fields to if dont exist
if ~isfield(meta,'qc'); meta.qc = ones(size(meta.f)); end;
if ~isfield(meta,'rmse'); meta.rmse=nan(size(meta.f));end
if ~isfield(meta,'dtrans'); meta.dtrans=nan(length(meta.f),3); end
if ~isfield(meta,'overlap'); meta.overlap=zeros(size(meta.f)); end % number existing data points overlapping file

while length(meta.f) >= 1
    
    fprintf('%d strips remaining\n',length(meta.f));
    
    % loop through files and record number of nan and non-nan mosaic grid
    % points overlap each strip.
    
    % only find overlapping data if any data exists in the tile
    meta.coregTestFlag=false(size(meta.f)); % flag which files have updated overlap to retest coregistration
    
    if any(N(:))
        
        % subset the DEMs to those that just overlap the rectangle of
        % existing data coverage
        anyN=any(N);
        minNx = x(find(anyN,1,'first'));
        maxNx = x(find(anyN,1,'last'));
        anyN=any(N');
        maxNy=y(find(anyN,1,'first'));
        minNy =y(find(anyN,1,'last'));
        clear anyN
        
        % find files within boundary
        n = find(meta.xmax > minNx & meta.xmin < maxNx & ...
            meta.ymax > minNy & meta.ymin < maxNy);
        
        % loop through files in boundary and find # of new/existing points
        overlap=meta.overlap;
        for i=1:length(n)
            
            % subset of N at data points in this file
            Nsub=N(meta.gridPointInd{n(i)});
            
            % count existing data points
            overlap(n(i))=sum(Nsub);
            
            % count new data points
            meta.gridPointN(n(i)) = sum(~Nsub);
            
        end
        
        % flag which files have updated overlap to retest coregistration
        meta.coregTestFlag = meta.overlap ~= overlap;
        
        % refresh overlap
        meta.overlap=overlap; overlap=[];
        
        % remove rendundant files (with already 100% coverage)
        redundantFlag = meta.gridPointN==0;
        
        if any(redundantFlag)
            
            % append redundant file records to ouput
            [m,meta] = appendRedundant(m,meta,redundantFlag);
            
            % remove redundant files from lists
            meta = structfun(@(x) ( x(~redundantFlag,:) ), meta, 'UniformOutput', false);
            
            fprintf('%d redundant files removed\n',sum(redundantFlag));
            
            if isempty(meta.overlap); break; end
        end
        
    end
    
    
    %% Test for coregistration rmse to use a selection criteria
    if ~disableCoregTest % check to see if disabled
        n=find(meta.coregTestFlag); %find data flagged for updating stats
        i= 1; % set iteration variable
        for i=1:length(n)
            
            fprintf('testing %d of %d\n',i,length(n))
            
            % check for mask fields and make mask cell if exists
            if isfield(meta,'maskPolyx')
                mask=[meta.maskPolyx{n(i)}(:) meta.maskPolyy{n(i)}(:)];
            else
                mask=cell(1,2);
            end
            
            [meta.dtrans(n(i),:),meta.rmse(n(i))] = coregister2tile(meta.f{n(i)},m,N,'mask',mask);
            
            if isnan(meta.rmse(n(i))); meta.overlap(n(i)) = 0; end
            
        end
        
    else
        meta.rmse(meta.overlap > 0)=0;
        meta.rmse(meta.overlap == 0)= NaN;
    end
    
    % while loop attempts to add data and will repeat if addition failure
    % results in meta data deletion, so that resorting of the data rank is
    % needed. Could probably be removed with a different indexing scheme.
    addFlag = false;
    while ~addFlag
        
        % buid selection array - need to invert meta.gridPointN account since ascending.
        A = [meta.rmse,double(meta.qc),max(meta.gridPointN)-meta.gridPointN,meta.overlap,(1:length(meta.f))'];
        
        % sort each row in priority order of columns
        A = sortrows(A);
        
        % while loop attempts to add data and will continue to add in A
        % list order as long as the meta index is not altered (e.g. in the
        % case where there is not enough coregistration overlap), so that
        % the data is just skipped but not removed from the meta list so
        % reording is not needed.
        skipFlag=false;
        j=1;
        while ~skipFlag && (j <= size(A,1))
            
            addErrors = true;
            
            % get the top selection index
            i = A(j,end);
            
            % if meta.rmse isnan, then this is a new cluster.
            if isnan(meta.rmse(i))
                meta.rmse(i)=0;
                meta.dtrans(i,:) = [0,0,0];
                c=c+1;
                addErrors=false;
            end
            
            if isfield(meta,'maskPolyx')
                mask=[meta.maskPolyx{i}(:) meta.maskPolyy{i}(:)];
            else
                mask=cell(1,2);
            end
            
            fprintf('adding strip: %s\n',meta.f{i})
            fprintf(...
                'quality: %d, coreg. RSME: %.2f, new data: %d points, co.cluster: %d\n',...
                meta.qc(i),meta.rmse(i),meta.gridPointN(i),c)
            
            [N,skipFlag] = addStrip2Mosaic( meta.f{i},m,meta.stripDate(i)-dy0,N,meta.dtrans(i,:)',meta.rmse(i),...
                'mergeMethod',mergeMethod,...
                'minNewPixels',minNewPixels,...
                'minOverlapPixels',minOverlapPixels,...
                'mask',mask,...
                'addErrors',addErrors);
            
            j=j+1;
        end
        
        % if skipFlag is still 0, none of the remaining DEMs will
        % coregister
        if ~skipFlag
            % clear meta.f to break outer loop
            meta.f=[];
            break;
        end
        
        % add this file name
        m.f=[m.f,meta.f(i)];
        
        % add meta data
        m.overlap   = [m.overlap,meta.overlap(i)];
        m.stripDate = [m.stripDate,meta.stripDate(i)];
        m.qcflag    = [m.qcflag,meta.qc(i)];
        m.maskPolyx = [m.maskPolyx,meta.maskPolyx(i)];
        m.maskPolyy = [m.maskPolyy,meta.maskPolyy(i)];
        
        if meta.rmse(i) == 0
            m.reg=[m.reg,{'N'}];
        else
            m.reg=[m.reg,{'A'}];
        end
        
        % remove this file from the meta struct
        in=1:length(meta.f); in(i)=[];
        meta = structfun(@(x) ( x(in,:) ), meta, 'UniformOutput', false);
        
        % check for success
        addFlag=~isnan(m.rmse(1,end));
        
    end
    % update cluster array
    C(C == 0 & N > 0)=c;
    m.C=C;
    
    % update data coverage field N in file
    m.N = N;
    
end

%% Generate Meta File
tileMeta(m);

function [m,meta] = appendRedundant(m,meta,redundantFlag)

% append redundant file records to ouput
m.f=[m.f,meta.f(redundantFlag)'];
m.reg=[m.reg,repmat({'-'},1,sum(redundantFlag))];
m.dtrans=[m.dtrans,nan(3,sum(redundantFlag))];
m.rmse=[m.rmse,nan(1,sum(redundantFlag))];
m.overlap = [m.overlap,meta.overlap(redundantFlag)'];
m.qcflag = [m.qcflag,meta.qc(redundantFlag)'];
m.maskPolyx =  [m.maskPolyx,meta.maskPolyx(redundantFlag)'];
m.maskPolyy =  [m.maskPolyy,meta.maskPolyy(redundantFlag)'];
m.stripDate = [m.stripDate,meta.stripDate(redundantFlag)'];
