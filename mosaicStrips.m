function mosaicStrips(meta,tiles,res,outdir,projstr,varargin)
% MOSAICSTRIPS upper level function for mosaicking DEM strips to a single tile, skips align, merge and tif export steps
%
%  mosaicStrips(db,tiles,res,outdir)
%  mosaicStrips(db,tiles,res,outdir,gcp)

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
        %strips2tile(meta,tiles.x0(n),tiles.x1(n),tiles.y0(n),tiles.y1(n),res,...
        %    outname,'disableReg','disableCoregTest');
        strips2tile(meta,tiles.x0(n),tiles.x1(n),tiles.y0(n),tiles.y1(n),res,...
            outname,'disableReg');
    end
    
    if exist(outname,'file')
        if exist(strrep(outname,'dem.mat','reg_dem.mat'),'file')
            fprintf('tile reg dem %s exists\n',tiles.I{n});
        else
            % Apply gcp registration if gcp's provided
            if exist('gcp','var')
                fprintf('applying ground control\n');
                m = matfile(outname);
                registerTile(m,gcp);
                clear m
            end
        end
    end
end

% tilef=dir([outdir,'/*m_dem.mat']);
% tilef = cellfun(@(x) [outdir,'/',x], {tilef.name}, 'UniformOutput',false);
% for i=1:length(tilef)
%     
%     fprintf('writing tif %d of %d\n',i,length(tilef))
%     
%     if exist(strrep(tilef{i},'dem.mat','reg_dem.mat'),'file')
%         fi=strrep(tilef{i},'dem.mat','reg_dem.mat');
%     else
%         fi=tilef{i};
%     end
%     
%     fprintf('source: %s\n',fi);
%     
%     load(fi,'x','y');
%     % crop buffer tile
%     x=x(101:end-100);
%     y=y(101:end-100);
%     
%     OutDemName = strrep(fi,'.mat','.tif');
%     if ~exist(OutDemName,'file');
%         
%         load(fi,'z');
%         z=z(101:end-100,101:end-100);
%         z(isnan(z)) = -9999;
%         
%         writeGeotiff(OutDemName,x,y,z,4,-9999,projstr)
%         clear z
%     end
%     
%     OutMatchtagName = strrep(fi,'dem.mat','matchtag.tif');
%     if ~exist(OutMatchtagName,'file');
%         load(fi,'mt');
%         mt =mt(101:end-100,101:end-100);
%         writeGeotiff(OutMatchtagName,x,y,mt,1,0,projstr)
%         clear mt
%     end
%     
% %     OutOrthoName = strrep(fi,'dem.mat','ortho.tif');
% %     if ~exist(OutOrthoName ,'file');
% %         load(fi,'or');
% %         or =or(101:end-100,101:end-100);
% %         writeGeotiff(OutOrthoName ,x,y,or,2,0,projstr)
% %         clear or
% %     end;
% %     
% %     OutDaynumName = strrep(fi,'dem.mat','daynum.tif');
% %     if ~exist(OutDaynumName,'file');
% %         load(fi,'dy');
% %         dy =dy(101:end-100,101:end-100);
% %         writeGeotiff(OutDaynumName,x,y,dy,2,0,projstr)
% %         clear dy
% %     end
% 
%     clear x y
% end
