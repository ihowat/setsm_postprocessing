function strips2tileRebuild(varargin)
% STRIPS2TILEREBUILD rebuild DEM mosaic tile from file
%
% strips2tileRebuild(m,res) where m is a mosaic matfile handle and res is
% the target resolution. A new file *resm_rebuild.mat will be saved to the
% same directory as m containing the new mosaic. If res is the same as the
% m resolution, the file will be identical aside from any changes made in
% the script. If the res is different, the new resolution will be written.
% Currently the translation statistics for m are used and the DEMs are not
% re-coregistered.

% Ian Howat, Ohio State University.

%% Set Parameters & Parse args
% set Y2K as day 0 for the daynumber grid
dy0 = datenum('jan 1 2000');

% parse inputs
m0=varargin{1};
res = varargin{2};

% get old grid spacing
res0 = m0.x(1,2)-m0.x(1,1);

% make outname
outname=strrep(m0.Properties.Source,'.mat','_rebuild.mat');
outname=strrep(outname,['_',num2str(res0),'m_'],['_',num2str(res),'m_']);

% get tile boundarys
tilex0 = min(m0.x);
tilex1 = max(m0.x);
tiley0 = min(m0.y);
tiley1 = max(m0.y);

% intialize outputs
x=[]; y=[]; z=[]; mt=[];  or=[];  dy=[];  f=[];  dtrans=[];  rmse=[];


%% Mosaic Grid Definition and Initialization
% define mosaic coorinate grid from tile boundaries
x = tilex0: res:tilex1;
y = tiley1:-res:tiley0;
y = y(:);

% make output file
save(outname,'x','y','-v7.3');
m = matfile(outname,'Writable',true);

tileVersion     = m0.version;
mergeMethod     = m0.mergeMethod;
disableReg      = m0.disableReg;
disableCoregTest= m0.disableCoregTest;
mergeMethodReg  = m0.mergeMethodReg;

m.version           = tileVersion;
m.mergeMethod       = mergeMethod;
m.disableReg        = disableReg;
m.disableCoregTest  = disableCoregTest;
m.mergeMethodReg    = mergeMethodReg;

% initialize mosaic grids
m.z = nan(length(y),length(x),'single'); % elevation data grid
m.ze = nan(length(y),length(x),'single'); % elevation error data grid
m.or =zeros(length(y),length(x),'int16'); % ortho imagery grid
m.mt = zeros(length(y),length(x),'uint8'); % matchtag data grid
m.dy = zeros(length(y),length(x),'int16'); % strip index grid
C = zeros(length(y),length(x),'int16'); % coregistration cluster grid
m.C=C;
N= zeros(length(y),length(x),'uint8'); % pixel data count
m.N=N;

% initialize mosaic meta data variables
m.f= []; m.reg=[]; m.overlap=[]; m.rmse=[]; m.dtrans=[];

% initialize cluster counter
c=1;

% extract meta data and copy to new file
f=m0.f;
m.f=f;
m.overlap = m0.overlap;
reg=m0.reg;
m.reg = reg;
m.stripDate = m0.stripDate;

% extract mask data if exists
a=whos(m0);
if any(strcmp({a.name},'maskPolyx'))
    maskPolyx=m0.maskPolyx;
    maskPolyy=m0.maskPolyy;
end

m.maskPolyx = maskPolyx;
m.maskPolyy = maskPolyy;

%% Add Registered Strips as Anchors
if ~disableReg && any(strcmp(reg,'R'))

    for i=1:length(reg)
        
        if ~strcmpi(reg{i},'R'); continue; end
        
        fprintf('adding anchor strip %d of %d: %s\n',i,length(reg),f{i})
        if exist('maskPolyx','var')
            mask=[maskPolyx{i}(:) maskPolyy{i}(:)];
        else
            mask=cell(1,2);
        end
        
        rmse = m0.rmse(1,i);
        dtrans = m0.dtrans(:,i);
        
        
        N = addStrip2Mosaic(f{i},m,m0.stripDate(1,i)-dy0,N,dtrans',rmse,...
            'mergeMethod',mergeMethodReg,...
            'mask',mask,...
            'addErrors',false);
        
        m.N = N;
        
        C(N ~= 0)=int16(1);
        m.C = C;
        m.Nreg = C;
        
    end
    
end

i=1;
for i=1:length(m0.rmse)
    
    addErrors=true;
    refineReg=false;

    if strcmpi(reg{i},'R'); continue; end
        
    if  isnan(m0.rmse(1,i))
        
        m.dtrans=[m.dtrans,[NaN;NaN;NaN]];
        m.rmse=[m.rmse,NaN];
        
        continue
    end

    if exist('maskPolyx','var')
        mask=[maskPolyx{i}(:) maskPolyy{i}(:)];
    else
        mask=cell(1,2);
    end

    % if m0.rmse and m0.dtrans is zero, this is a new cluster
    if m0.rmse(1,i) == 0 && ~any(m0.dtrans(:,i) ~= 0)
        c=c+1;
        addErrors=false;
        dtrans=[0;0;0];
        rmse=0;
    else
        rmse = m0.rmse(1,i);
        dtrans = m0.dtrans(:,i);
        
        % refine coregistration if different size grid. Note it maye be
        % worthwhile to just turn this off to save time as Im not sure this
        %really improves the registration that much. 
        if res ~= res0; refineReg = true; end;
        
    end
    
    fprintf('adding strip: %s\n',f{i})
    fprintf(...
        'quality: %d, coreg. RSME: %.2f, co.cluster: %d\n',...
        m0.qcflag(1,i),rmse,c)
    
    [N,~] = addStrip2Mosaic(f{i},m,m0.stripDate(1,i)-dy0,N,dtrans',rmse,...
        'mergeMethod',mergeMethod,...
        'mask',mask,...
        'addErrors',addErrors,...
        'refineReg',refineReg);
    
    % update cluster array
    C(C == 0 & N > 0)=c;
    m.C=C;

    % update data coverage field N in file
    m.N = N;
    %
end
%
% %% Generate Meta File
tileMeta(m);

