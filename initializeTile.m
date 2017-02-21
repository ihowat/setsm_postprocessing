function [x,y,c,C,N,m,meta] = initializeTile(meta,tilex0,tilex1,tiley0,tiley1,res,outname,tileVersion,disableReg, disableCoregTest,mergeMethod,mergeMethodReg);

% initializeTile mosaic tile grid definition and initialization or restart
%
% [x,y,c,C,N,m] = initializeTile(...
%       meta,tilex0,tilex1,tiley0,tiley1,res,outname) 
%       ARGINS: meta is the metastruct, tile* are the tile grid ranges,
%       res is the output grid resolution and outname is the output file 
%       name. 
%       ARGOOUTS: x,y are grid coordinate vectors, C is the coregistration
%       cluster number array, c is the initial coregisteration cluster
%       number, N is the grid point data flag and m is the output file
%       handle.
%
%   Ian Howat, ihowat@gmail.com, Ohio State

% Initialize new mosaic file if not exists
if ~exist(outname,'file')
    
    % define mosaic coorinate grid. Add a 100-pixel buffer for aligning/merging
    % the tiles
    x = tilex0-100*res: res:tilex1+100*res;
    y = tiley1+100*res:-res:tiley0-100*res;
    y = y(:);
    
    % make output file
    save(outname,'x','y','-v7.3');
    m = matfile(outname,'Writable',true);
    
    m.version=tileVersion;
    m.mergeMethod=mergeMethod;
    m.disableReg=disableReg;
    m.disableCoregTest=disableCoregTest;
    m.mergeMethodReg=mergeMethodReg;

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
    m.qcflag =  []; m.maskPolyx =  []; m.maskPolyy =  [];  m.stripDate = [];
    
    % initialize cluster counter
    c=1;

% Attempt to restart if output file exists
else %file already exists so try and pick up where it left off
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