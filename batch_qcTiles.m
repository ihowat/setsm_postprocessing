function batch_qcTiles(tileNames)

FIGURE_CLEANUP = onCleanup(@() eval('close all'));

%tileFile = '/data4/REMA/region_06_luitpold_coast/mosaic_reg_qc_feather2/40m/36_18_40m_dem.mat';
dbasefile  = 'V:\pgc\data\scratch\claire\repos\setsm_postprocessing_pgc\arcticDEMdatabase3_2m_unf_20200519.mat'; % database file
tileDir = 'V:\pgc\data\elev\dem\setsm\ArcticDEM\mosaic\2m_v3.2_tileqc'; %directory containing tiles
%tileName = '31_39';

dbasedir_local = [getenv('USERPROFILE'),'\setsm_postprocessing_dbase'];
[dbasefile] = copy_dbase_local(dbasedir_local, dbasefile);

tileNames = strsplit(tileNames,',');

for i=1:length(tileNames)
    tileFiles{i} = [tileDir,'\',char(tileNames(i)),'\',char(tileNames(i)),'_8m_dem.mat'];
end
%fileDates = [tileFiles.datenum];
%tileFiles = {tileFiles.name};
%tileFiles = cellfun( @(x) [tileDir,'\',tileName,'\',x], tileFiles,'uniformoutput',0);

% use this if to only look at recently created files
% n = fileDates > datenum('31-Jan-2018 12:00:00'); 
% tileFiles = tileFiles(n);

fprintf('Loading db\n');
meta=load(dbasefile);
f1=figure('Name','BQCT');

for i=1:length(tileFiles)

    if exist(tileFiles{i},'file')
        fprintf('tile %s, %d of %d\n',tileFiles{i},i,length(tileFiles));
        [~,tn,~] = fileparts(tileFiles{i});
        tn=strrep(tn,'_8m_dem','');
        f1.Name = tn;
        qcTile(tileFiles{i},meta);
    end

end
end


function [dstfile] = make_local_copy(srcfile, dstdir)
    srcfile_stats = dir(srcfile);
    srcfname = srcfile_stats.name;
    dstfile = [dstdir,'\',srcfname];
    
    if exist(dstdir, 'dir') ~= 7
        fprintf('Making local copy of db file at %s ...', dstfile);
        mkdir(dstdir);
        
    elseif exist(dstfile, 'file') == 2
        dstfile_stats = dir(dstfile);
        if dstfile_stats.datenum == srcfile_stats.datenum
            % Local copy doesn't needed to be updated.
            return;
        else
            fprintf('Updating local copy of db file at %s ...', dstfile);
        end
        
    else
        fprintf('Making local copy of db file at %s ...', dstfile);
    end
    
    status = copyfile(srcfile, dstdir);
    if status == 0
        fprintf(' failed!\n');
        dstfile = '';
    elseif status == 1
        fprintf(' success!\n');
    end
end


function [varargout] = copy_dbase_local(dbasedir_local, varargin)
    for i = 1:length(varargin)
        dbasefile = varargin{i};
        dbasefile_local = make_local_copy(dbasefile, dbasedir_local);
        if isempty(dbasefile_local)
            fprintf('Falling back to network load\n');
        else
            dbasefile = dbasefile_local;
        end
        varargout{i} = dbasefile;
    end
end
