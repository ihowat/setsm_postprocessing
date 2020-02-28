function generate_qc_files_and_setfacl(regionRootDir)

region_dir_prefix = 'arcticdem_';
region_subdir = 'strips_unf/2m';

d = dir(regionRootDir);
regionDirs = d(cellfun(@(x) startsWith(x, region_dir_prefix), {d.name}));
if isempty(regionDirs)
    fprintf("No region folders found in root dir matching prefix '%s'\n", region_dir_prefix);
    return;
end
regionDirs = {regionDirs.name};

for i = 1:numel(regionDirs)

    regionFolder = char(regionDirs(i));
    
    stripDir = [regionRootDir,'/',regionFolder,'/',region_subdir];
    if exist(stripDir,'dir') ~= 7
        fprintf("%s does not have a subdir '%s', skipping\n", regionFolder, region_subdir);
        continue;
    end
    
    qcFile = [stripDir,'/','qc.mat'];
    qcLockFile = strrep(qcFile, 'qc.mat', 'qc.mat.lock');
    
    
    run_setfacl = true;
    if exist(qcFile, 'file')
        
        cmd = sprintf('getfacl %s', qcFile);
        [status, out] = system(cmd);
        if status ~= 0
            error(out);
        end
        if contains(out, 'group:esci-pgc-users:rwx')
            run_setfacl = false;
        end
        
    else
    
        fprintf('Creating %s\n', qcFile);
        
        fileNamesTemp = dir([stripDir,'/*/*dem_10m_shade_masked.tif']);
        fileNames = fullfile({fileNamesTemp.folder},{fileNamesTemp.name});
        %fileNames = cellfun(@(x) [x.folder,'/',x.name], {fileNamesTemp}, 'uniformOutput',false);
        fileNames = fileNames(:);

        flag = zeros(size(fileNames));
        x = cell(size(fileNames));
        y = cell(size(fileNames));

        save(qcFile, 'fileNames','x','y','flag');
        
    end
    
    
    if run_setfacl
        cmd = sprintf('setfacl -m g:esci-pgc-users:rwx %s', qcFile);
        [status, out] = system(cmd);
        if status ~= 0
            error(out);
        end
    end
    
    
    run_setfacl = true;
    if exist(qcLockFile, 'file')
        
        cmd = sprintf('getfacl %s', qcLockFile);
        [status, out] = system(cmd);
        if status ~= 0
            error(out);
        end
        if contains(out, 'group:esci-pgc-users:rwx')
            run_setfacl = false;
        end
        
    else
        
        fprintf('Creating %s\n', qcLockFile);
        
        cmd = sprintf('touch %s', qcLockFile);
        [status, out] = system(cmd);
        if status ~= 0
            error(out);
        end
        
    end
    
    if run_setfacl
        cmd = sprintf('setfacl -m g:esci-pgc-users:rwx %s', qcLockFile);
        [status, out] = system(cmd);
        if status ~= 0
            error(out);
        end
    end
    
end

end
