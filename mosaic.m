function mosaic(db,tiles,res,outdir,varargin)
% MOSAIC upper level function for mosaicking DEM strips to a set of tiles
%
%  mosaic(db,tiles,res,outdir)
%  mosaic(db,tiles,res,outdir,gcp)

%% parse and error check varargin
if ~isempty(varargin)
       
       gcp=varargin{1};
       
       if ~isfield(gcp,'x') || ~isfield(gcp,'y') || ~isfield(gcp,'z')
           error('gcp arg must be structure with fields x,y,z')
       end
end

tilef = cellstr([char(tiles.I),...
    repmat(['_',num2str(res),'m_dem.mat'],length(tiles.I),1)]);

%% Loop through tiles and mosaic
for n=1:length(tiles.I)
    
    outname=[outdir,'/',tiles.I{n},'_',num2str(res),'m_dem.mat'];
    
    if exist(outname,'file')
        fprintf('tile dem %s exists\n',tiles.I{n})
    %                 fprintf('tile %s exists deleting\n',tiles.I{n})
    %                 delete(outname);
    else
        fprintf('building tile %s: %d of %d\n',tiles.I{n},n,length(tiles.I))
        strips2tile(db,tiles.x0(n),tiles.x1(n),tiles.y0(n),tiles.y1(n),res,...
            outname,'disableReg');
    end
    
    if exist(outname,'file')
        if exist(strrep(outname,'dem.mat','reg_dem.mat'),'file')
            fprintf('tile reg dem %s exists\n',tiles.I{n});
        else
            % Apply gcp registration if gcp's provided
            if exist('gcp','var')
                fprintf('appling ground control\n');
                m = matfile(outname);
                registerTile(m,gcp);
                clear m
            end
        end
    end
end


%% Align unregistered tiles

% find unregistered tiles by intersect/removing with registered list
f=dir([outdir,'/*m_dem.mat']);
f=cellfun(@(x) [outdir,'/',x], {f.name}, 'UniformOutput',false);

freg=dir([outdir,'/*m_dem_reg.mat']);
freg=cellfun(@(x) [outdir,'/',x], {freg.name}, 'UniformOutput',false);

freg=strrep(freg,'_reg','');
[~,IA]=intersect(f,freg);
f(IA)=[];

batchAlignTile(f,freg);
 

%% Blend tile edges
% get list of tiles from current directory
tilef=dir([outdir,'/*m_dem_reg.mat']);
tilef = cellfun(@(x) [outdir,'/',x], {tilef.name}, 'UniformOutput',false);

batchMergeTileBuffer(tilef);

%% Crop Buffers and Write Tiles To Geotiffs
for i=1:length(tilef)
    
    fprintf('writing tif %d of %d\n',i,length(tilef))
    
    if exist(strrep(tilef{i},'dem.mat','reg_dem.mat'),'file')
        fi=strrep(tilef{i},'dem.mat','reg_dem.mat');
    else
        fi=tilef{i};
    end
    
    fprintf('source: %s\n',fi);
    
    load(fi,'x','y');
    % crop buffer tile
    x=x(101:end-100);
    y=y(101:end-100);
    
    OutDemName = strrep(fi,'.mat','.tif');
    if ~exist(OutDemName,'file');
        
        load(fi,'z');
        z=z(101:end-100,101:end-100);
        z(isnan(z)) = -9999;
        
        writeGeotiff(OutDemName,x,y,z,4,-9999,'polar stereo north')
        clear z
    end
    
    OutMatchtagName = strrep(fi,'dem.mat','matchtag.tif');
    if ~exist(OutMatchtagName,'file');
        load(fi,'mt');
        mt =mt(101:end-100,101:end-100);
        writeGeotiff(OutMatchtagName,x,y,mt,1,0,'polar stereo north')
        clear mt
    end
    
    OutOrthoName = strrep(fi,'dem.mat','ortho.tif');
    if ~exist(OutOrthoName ,'file');
        load(fi,'or');
        or =or(101:end-100,101:end-100);
        writeGeotiff(OutOrthoName ,x,y,or,2,0,'polar stereo north')
        clear or
    end;
    
    OutDaynumName = strrep(fi,'dem.mat','daynum.tif');
    if ~exist(OutDaynumName,'file');
        load(fi,'dy');
        dy =dy(101:end-100,101:end-100);
        writeGeotiff(OutDaynumName,x,y,dy,2,0,'polar stereo north')
        clear dy
    end

    clear x y
end
