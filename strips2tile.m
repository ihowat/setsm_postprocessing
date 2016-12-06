function strips2tile(meta,tilex0, tilex1, tiley0, tiley1,res,outname,varargin)
% STRIPS2TILE mosaic strips to rectangular tile
%
%[x,y,z,mt,or,dy,f,dtrans,rmse] = ...
%           strips2tile(meta,tilex0, tilex1, tiley0, tiley1,res,outname)
%   
% where meta is the database structure, tilex and tiley are the coordinate
% ranges of the tile, res is the tile resolution in m and outname is the
% tile output file name.
%
%[...] = strips2tile(...,'options') where the following options are
%availble:
%       'disableReg' disables use of any existing GCP registration
%       files/data
%       'disableCoregTest' disables the coregistraiton optimization test
%       to save time
%       'mergeMethod','methodstring', specifies the mergeing method as
%       'feather' (default), 'underprint' or 'warp'.
%
%   subfunctions: readreg, addStrip2Mosaic 
%
%   Ian Howat, Ohio State University
%   Version 1.0; 28-Jul-2016 09:50:32 (beta versions preceeded this)
%   Version 2.0; 8-Oct-2016 15:19:52
%   - preferentially stacks in order of best coregistration fit (if
%     disabled, sorts by quality rank and then new data coverage)
%	- added redundant data check & skip
%   - faster existing data search using gridPointInd field
%   - allows for use of a qc rank list
%   Version 2.1; 06-Dec-2016 14:22:22
%   - fixed bug in which disableCoregTest resuted in no coregistration in
%   addStrip2Mosiac

tileVersion='2.1';

%% Set Parameters & Parse args
% set Y2K as day 0 for the daynumber grid
dy0 = datenum('jan 1 2000');

% option to disable using control registration files
disableReg=false;

% option to disable coregistration optimization test
disableCoregTest=false;

% default merging method
mergeMethod='feather';

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

end

fprintf('Using merge method: %s\n',mergeMethod)

%% Initialize new mosaick file if not exists
if ~exist(outname,'file')
    
    % intialize outputs
    x=[]; y=[]; z=[]; mt=[];  or=[];  dy=[];  f=[];  dtrans=[];  rmse=[];
    
    % make tile boundary polygon
    tilevx = [tilex0;tilex0;tilex1;tilex1;tilex0];
    tilevy = [tiley0;tiley1;tiley1;tiley0;tiley0];
    
    %% Search for strips within this tile
    
    % quick search: find strips within range of this tile. This does not
    % account for background area of around strips but just pairs them down to
    % speed the poly intersection loop
    n = meta.xmax > tilex0 & meta.xmin < tilex1 & ...
        meta.ymax > tiley0 & meta.ymin < tiley1;
    
    % if no overlap, return
    if ~any(n); fprintf('no strip overlap\n'); return; end
    
    % par down database structure to possible overlapping tiles
    meta = structfun(@(x) ( x(n) ), meta, 'UniformOutput', false);
    
    % search for all strip footprints overlapping this tile
    n=zeros(size(meta.f));
    for i=1:length(n)
        n(i) = any(inpolygon(meta.x{i},meta.y{i},tilevx,tilevy)) | ...
            any(inpolygon(tilevx,tilevy,meta.x{i},meta.y{i}));
    end
    
    % if no overlap, return
    if ~any(n); fprintf('no strip overlap'); return; end
    
    % remove files with no overlap
    meta = structfun(@(x) ( x(logical(n)) ), meta, 'UniformOutput', false);
    
    %% Mosaic Grid Definition and Initialization
    
    % define mosaic coorinate grid. Add a 100-pixel buffer for aligning/merging
    % the tiles
    x = tilex0-100*res: res:tilex1+100*res;
    y = tiley1+100*res:-res:tiley0-100*res;
    y = y(:);
    
    % make output file
    save(outname,'x','y','-v7.3');
    m = matfile(outname,'Writable',true);
    
    % initialize mosaic grids
    m.z = nan(length(y),length(x),'single'); % elevation data grid
    m.or =zeros(length(y),length(x),'int16'); % ortho imagery grid
    m.mt = zeros(length(y),length(x),'uint8'); % matchtag data grid
    m.dy = zeros(length(y),length(x),'int16'); % strip index grid
    C = zeros(length(y),length(x),'int16'); % coregistration cluster grid
    m.C=C;
    N= zeros(length(y),length(x),'uint8'); % pixel data count
    m.N=N;
    
    % initialize mosaic meta data variables
    m.f= []; m.overlap=[]; m.rmse=[]; m.dtrans=[];
    
    % initialize cluster counter
    c=1;

%% Attempt to restart if output file exists
else %file already exists so try and pick up where it left off
    % WARNING not tested - may not pick up where you want it to
    fprintf('reading existing file and trying to pick up where it left off\n')
    m = matfile(outname,'Writable',true);
    
    if length(m.rmse) - length(m.f) == 1
        rmse=m.rmse; rmse(end) = []; m.rmse=rmse; rmse = [];
        dtrans=m.dtrans; dtrans(:,end) = []; m.dtrans=dtrans; dtrans = [];
    elseif  length(m.rmse) - length(m.f) > 1
         error('variable sizes in the existing mosaic file not reconcileable')
    end
     
    % remove already added files
    [~,IA]= intersect(meta.f,m.f);
    in=1:length(meta.f); in(IA)=[];
    meta = structfun(@(x) ( x(in,:) ), meta, 'UniformOutput', false);
    
    x=m.x;
    y=m.y;
    N=m.N;
    C=m.C;
    c=max(C(:));  
end

%% Compile Registration Data
if ~disableReg % check for disable registration argument
    regfiles=strrep(meta.f,'meta.txt','reg.txt');
    nreg= find( cellfun( @exist, regfiles) == 2);
    
    meta.dtrans=nan(size(meta.f,1),3);
    meta.rmse=nan(size(meta.f));
    
    if ~isempty(nreg)
        regfiles=regfiles(nreg);
        
        [~,~,~,dtrans,pctile]=cellfun( @readreg,regfiles,'uniformoutput',0);
        pctile=cell2mat(pctile);
        dtrans=cell2mat(dtrans);
        p75=pctile(:,6);
        
        meta.dtrans(nreg,:)=dtrans;
        meta.rmse(nreg)=p75;
    end
    
    %filter out data over 4m std
    meta.dtrans(meta.rmse > 4,:) = NaN;
    meta.rmse(meta.rmse > 4) = NaN;
    
    % sort database by p75
    [~,n] = sort(meta.rmse,'ascend');
    meta = structfun(@(x) ( x(n,:) ), meta, 'UniformOutput', false);
    
    nreg=find(~isnan(meta.rmse));
    
    %% Add Registered Strips as Anchors
    if ~isempty(nreg)
        
        for i=1:length(nreg)
            
            m.overlap=[m.overlap,0];
            
            fprintf('adding anchor strip %d of %d: %s\n',i,length(nreg),meta.f{i})
            if isfield(meta,'maskPolyx')
                mask=[meta.maskPolyx{i}(:) meta.maskPolyy{i}(:)];
            else
                mask=cell(1,2);
            end
            N = addStrip2Mosaic( meta.f{i},m,meta.stripDate(i)-dy0,N,...
                meta.dtrans(i,:)',meta.rmse(i),...
                'mergeMethod','underprint',...
                'mask',mask);

            % add this file name
            m.f=[m.f,meta.f(i)];
            m.N = N;
            
        end

        % remove these files from the meta struct
        in=length(nreg)+1:length(meta.f);
        meta = structfun(@(x) ( x(in,:) ), meta, 'UniformOutput', false);
        
        C(N ~= 0)=int16(1);
        m.C = C;
        
    end

end

%% Add dummy qc field if doesnt exist
if ~isfield(meta,'qc'); meta.qc = ones(size(meta.f)); end;

%% Build Grid Point Index Field
% make a cell for each file containing the col-wise indices of the data
% points on the tile grid. This is used search for new or existing data on
% the grid within each search iteration.

% first make sure does not exist in case we're saving this info
if ~isfield(meta,'gridPointInd')
    
    fprintf('calculating tile pixel coverage for each strip\n');
    
    % initialize output field  
    meta.gridPointInd=cell(size(meta.f));
    
    % can be slow, so we'll use a counter
    count = 0;
    
    % file loop
    for i=1:length(meta.f)
        
        % counter
        if i>1 for p=1:count fprintf('\b'); end; %delete line before
            count = fprintf('strip %d of %d',i,size(meta.f,1));
        end
        
        % locate grid pixels within footprint polygon
        BW = roipoly(x, y, N, meta.x{i}, meta.y{i});
        
        % if mask exists, apply it
        if isfield(meta,'maskPolyx')
            mask=[meta.maskPolyx{i}(:) meta.maskPolyy{i}(:)];
            
            for j=1:size(mask,1)
                if ~isempty(mask{j,1}) && ~isempty(mask{j,2})
                    BW(roipoly(x,y,N,mask{j,1},mask{j,2}))=0;
                end
            end
            
        end
   
        % convert BW mask to col-wise indices and save to cell
        meta.gridPointInd{i}=find(BW);
        
        % get rid of mask
        clear BW
    end
    
    % clear counter line
    for p=1:count fprintf('\b'); end;
    
    % make another field with the number of data points
    meta.gridPointN=cellfun(@length,meta.gridPointInd);

end

%% Sequential Floating Strip Addition Loop
% Sequentially add strips to the mosaic by selecting the file with the most
% new data coverage, with enough overlap to coregister to exisiting data.

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
            m.f=[m.f,meta.f(redundantFlag)'];
            m.dtrans=[m.dtrans,nan(3,sum(redundantFlag))];
            m.rmse=[m.rmse,nan(1,sum(redundantFlag))];
            m.overlap = [m.overlap,meta.overlap(redundantFlag)'];
            
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
            
            % get the top selection index
            i = A(j,end);
            
            % if meta.rmse isnan, then this is a new cluster.
            if isnan(meta.rmse(i)) 
                meta.rmse(i)=0;
                meta.dtrans(i,:) = [0,0,0];
                c=c+1; 
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
                'mask',mask);
            
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
        
        % add overlap
        m.overlap = [m.overlap,meta.overlap(i)];
        
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
m.version=tileVersion;
tileMeta(m);
