function strips2tile(db,tilex0, tilex1, tiley0, tiley1,res,outname)
% STRIPS2TILE mosaic strips to rectangular tile
%
%[x,y,z,mt,or,dy,f,dtrans,rmse] = ...
%           strips2tile(db,tilex0, tilex1, tiley0, tiley1,res,includeList)
%   
% where db is the database structure, tilex and tiley are the coordinate
% ranges of the tile, res is the tile resolution in m and includeList is a
% list of strips to include (optional).

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
% set Y2K as day 0 for the daynumber grid
dy0 = datenum('jan 1 2000');

% initialize mosaic meta data variables
m.f= []; m.overlap=[]; m.rmse=[]; m.dtrans=[];

% initialize cluster counter
c=0;

% %% determine if reg.txt files exist
% regfiles=strrep(db.f,'meta.txt','reg.txt');
% nreg = find ( cellfun( @exist, regfiles) == 2);
% regfiles=regfiles(nreg);
% 
% [~,~,~,trans,pctile]=cellfun( @readreg,regfiles,'uniformoutput',0);
% pctile=cell2mat(pctile);
% trans=cell2mat(trans);
% p75=pctile(:,6);
% 
% %sort by 75th percentile of residuals
% [~,n]=sort(pctile(:,6),'ascend');
% 
% nreg = nreg(n);
% trans=trans(n,:);
% 
% i=1;
% for i=1:length(nreg);
%          
%     
%     % N = floatStrip2Mosaic_v3(db.f{i}, m, db.stripDate(i)-dy0,N);
%         
%      fprintf('adding registered strip: %s\n',db.f{nreg(i)})
%    
%      N = floatStrip2Mosaic_v3( db.f{nreg(i)}, m,...
%                             db.stripDate(nreg(i))-dy0, N,...
%                             trans(i,:));
%                         
%       m.N = N;
% end
%         
%  % remove these files from the meta struct
% in=1:length(db.f); in(nreg)=[];
%  db = structfun(@(x) ( x(in) ), db, 'UniformOutput', false);


%% Sequential Strip Addition Loop
% sequentially add strips to the mosaic. Currently adds by the most overlap
% with the tile

% calculate tile coverage for all strips
fprintf('calculating tile pixel coverage for each strip\n');
db.tileCoverage=zeros(size(db.f));
for i=1:length(db.f)
    BW = roipoly(x, y, N, db.x{i}, db.y{i});
    db.tileCoverage(i)=sum(BW(:));
    clear BW
end

while length(db.f) >= 1;
    
    fprintf('%d strips remaining\n',length(db.f));
    
    % loop through files and record number of nan and non-nan mosaic grid
    % points overlap each strip.
    Ndata=zeros(size(db.f));
    if any(N(:))
        anyN=any(N);
        minNx = m.x(1,find(anyN,1,'first')); 
        maxNx = m.x(1,find(anyN,1,'last'));
        anyN=any(N');
        maxNy=m.y(find(anyN,1,'first'),1); 
        minNy =m.y(find(anyN,1,'last'),1);
        clear anyN
        
        n = find(db.xmax > minNx & db.xmin < maxNx & ...
                db.ymax > minNy & db.ymin < maxNy);

        for i=1:length(n);
            BW = roipoly(x, y, N, db.x{n(i)}, db.y{n(i)});
            Ndata(n(i))=sum(N(BW(:))~=0);
            clear BW
        end
    end
    
    % if non-nan overlap exists select strip with most overlap to non-nans
    if any(Ndata)
        [m.overlap(1,length(m.f)+1),i]=max(Ndata);
        fprintf('adding floating strip: %s\n',db.f{i})
       
%         [z,mt,or,dy,dtrans(:,length(f)),rmse(length(f))] = floatStrip2Mosaic(...
%             db.f{i},x,y,z,mt,or,dy, db.stripDate(i)-dy0);

        N = floatStrip2Mosaic(db.f{i}, m, db.stripDate(i)-dy0,N);
        
    else
        % otherwise select strip with most mosaic overlap
        [~,i]=max(db.tileCoverage);
        m.overlap=[m.overlap,0];
        
        fprintf('adding anchor strip: %s\n',db.f{i})

        % start new coreg cluster
        c=c+1;
        
%         [z,mt,or,dy] = refStrip2Mosaic(...
%             db.f{i},x,y,z,mt,or,dy,db.stripDate(i)-dy0);
%         dtrans= [dtrans,zeros(3,1)];
%         rmse  = [rmse,0];    
        
          N = refStrip2Mosaic( db.f{i},m,db.stripDate(i)-dy0,N);
    end
    
    % add this file name 
    m.f=[m.f,db.f(i)];
    
    C(C == 0 & N > 0)=c;
    
    m.C=C;
    m.N = N;
    
    if length(db.f) == 1; break; end
    
    % remove this file from the meta struct
    in=1:length(db.f); in(i)=[];
    db = structfun(@(x) ( x(in) ), db, 'UniformOutput', false);
    
end

