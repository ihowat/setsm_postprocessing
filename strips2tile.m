function strips2tile(db,tilex0, tilex1, tiley0, tiley1,res,outname,varargin)
% STRIPS2TILE mosaic strips to rectangular tile
%
%[x,y,z,mt,or,dy,f,dtrans,rmse] = ...
%           strips2tile(db,tilex0, tilex1, tiley0, tiley1,res,outname)
%   
% where db is the database structure, tilex and tiley are the coordinate
% ranges of the tile, res is the tile resolution in m and outname is the
% tile output file name.
%
%   subfunctions: readreg, addStrip2Mosaic 
%
%   Ian Howat, Ohio State University
%   version 1; 28-Jul-2016 09:50:32 (beta versions preceeded this)

%  version number
tileVersion='1.0';

% set Y2K as day 0 for the daynumber grid
dy0 = datenum('jan 1 2000');

% option to disable using control registration files
disableReg=false;
if ~isempty(varargin); disableReg=true; end

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
    n = db.xmax > tilex0 & db.xmin < tilex1 & ...
    db.ymax > tiley0 & db.ymin < tiley1;

    % if no overlap, return
    if ~any(n); fprintf('no strip overlap\n'); return; end

    % par down database structure to possible overlapping tiles
    db = structfun(@(x) ( x(n) ), db, 'UniformOutput', false);

    % search for all strip footprints overlapping this tile
    n=zeros(size(db.f));
    for i=1:length(n);
        n(i) = any(inpolygon(db.x{i},db.y{i},tilevx,tilevy)) | ...
           any(inpolygon(tilevx,tilevy,db.x{i},db.y{i}));
    end

    % if no overlap, return
    if ~any(n); fprintf('no strip overlap'); return; end

    % remove files with no overlap
    db = structfun(@(x) ( x(logical(n)) ), db, 'UniformOutput', false);

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

else %file already exists so try and pick up where it left off
    % WARNING not tested - may not pick up where you want it to
    warning('reading existing file and trying to pick up where it left off, but not tested or recommended')
    m = matfile(outname,'Writable',true);
    [~,IA,IB]= intersect(db.f,m.f);
    in=1:length(db.f); in(IA)=[];
    db = structfun(@(x) ( x(in) ), db, 'UniformOutput', false);
  
    x=m.x;
    y=m.y;
    N=m.N;
    C=m.C;
    c=max(C(:));
end

%% Compile Registration Data
if ~disableReg % check for disable registration argument
    regfiles=strrep(db.f,'meta.txt','isreg.txt');
    nreg= find( cellfun( @exist, regfiles) == 2);

    db.trans=nan(size(db.f,1),3);
    db.p75=nan(size(db.f));

    if ~isempty(nreg)
        regfiles=regfiles(nreg);
  
        [~,~,~,trans,pctile]=cellfun( @readisreg,regfiles,'uniformoutput',0);
        pctile=cell2mat(pctile);
        trans=cell2mat(trans);
        p75=pctile(:,6);
    
        db.trans(nreg,:)=trans;
        db.p75(nreg)=p75;
    end

    %filter out data over 4m std
    db.p75(db.p75 > 4) = NaN;

    % sort database by p75
    [~,n] = sort(db.p75,'ascend');
    db = structfun(@(x) ( x(n,:) ), db, 'UniformOutput', false);

    nreg=find(~isnan(db.p75));

    %% Add Registered Strips as Anchors
    for i=1:length(nreg);
    
        m.overlap=[m.overlap,0];
        
        fprintf('adding anchor strip %d of %d: %s\n',i,length(nreg),db.f{i})
        
        N = addStrip2Mosaic( db.f{i},m,db.stripDate(i)-dy0,N,db.trans(i,:),'warp');
        
        % add this file name 
        m.f=[m.f,db.f(i)];
        
        m.N = N;     
    end

    % remove these files from the meta struct
    in=length(nreg)+1:length(db.f);
    db = structfun(@(x) ( x(in,:) ), db, 'UniformOutput', false);
    
end

%% Sequential Floating Strip Addition Loop
% sequentially add strips to the mosaic. Currently adds by the most overlap
% with the tile

while length(db.f) >= 1;
    
    fprintf('%d strips remaining\n',length(db.f));
    
    % loop through files and record number of nan and non-nan mosaic grid
    % points overlap each strip.
    Ndata=zeros(size(db.f));
    if any(N(:))
        anyN=any(N);
        minNx = x(find(anyN,1,'first'));
        maxNx = x(find(anyN,1,'last'));
        anyN=any(N');
        maxNy=y(find(anyN,1,'first'));
        minNy =y(find(anyN,1,'last'));
        clear anyN
        
        n = find(db.xmax > minNx & db.xmin < maxNx & ...
            db.ymax > minNy & db.ymin < maxNy);
        
        fprintf('calculating tile strip data overlap\n');
        count=0;
        for i=1:length(n);
            if i>1 for p=1:count fprintf('\b'); end; %delete line before
                count = fprintf('strip %d of %d',i,length(n));
            end
            
            BW = roipoly(x, y, N, db.x{n(i)}, db.y{n(i)});
            Ndata(n(i))=sum(N(BW(:))~=0);
            clear BW
        end
        for p=1:count fprintf('\b'); end;
    end
    
    % if non-nan overlap exists select strip with most overlap to non-nans
    if any(Ndata)
        [m.overlap(1,length(m.f)+1),i]=max(Ndata);
        fprintf('adding floating strip: %s\n',db.f{i})
        
        N = addStrip2Mosaic( db.f{i},m,db.stripDate(i)-dy0,N,[],'feather');
        
    else
         % otherwise select strip with most mosaic overlap
        if ~isfield('db','tileCoverage');
            % calculate tile coverage for all strips
            fprintf('calculating tile pixel coverage for each strip\n');
            db.tileCoverage=zeros(size(db.f));
            
            count = 0;
            for i=1:length(db.f)
               
                if i>1 for p=1:count fprintf('\b'); end; %delete line before
                    count = fprintf('strip %d of %d',i,size(db.f,1));
                end
                BW = roipoly(x, y, N, db.x{i}, db.y{i});
                db.tileCoverage(i)=sum(BW(:));
                clear BW
            end
            for p=1:count fprintf('\b'); end;
        end

        [~,i]=max(db.tileCoverage);
        m.overlap=[m.overlap,0];
        
        fprintf('adding anchor strip: %s\n',db.f{i})
        
        % start new coreg cluster
        c=c+1;
        
       N = addStrip2Mosaic( db.f{i},m,db.stripDate(i)-dy0,N,[0,0,0],'feather');
       
    end

    % add this file name
    m.f=[m.f,db.f(i)];
    
    C(C == 0 & N > 0)=c;
    
    m.C=C;
    m.N = N;
    
    if length(db.f) == 1; break; end
    
    % remove this file from the meta struct
    in=1:length(db.f); in(i)=[];
    db = structfun(@(x) ( x(in,:) ), db, 'UniformOutput', false);
    
end

%% Generate Meta File
m.version=tileVersion;
tileMeta(m);

