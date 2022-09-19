function generate_qc_files_and_setfacl(regionRootDir)

region_dir_prefix = 'arcticdem_';
region_subdir = 'strips_unf/2m';

%%% CHECK THIS SETTING %%%
report_number_of_strips_to_append_but_dont_actually_append = true;
%%% CHECK THIS SETTING %%%

d = dir(regionRootDir);
regionDirs = d(cellfun(@(x) startsWith(x, region_dir_prefix), {d.name}));
if isempty(regionDirs)
    fprintf('No region folders found in root dir matching prefix "%s"\n', region_dir_prefix);
    return;
end
regionDirs = {regionDirs.name};

for i = 1:numel(regionDirs)

    regionFolder = char(regionDirs(i));
%    if ~strcmp(regionFolder, 'arcticdem_05_greenland_northeast'); continue; end
%    if ~strcmp(regionFolder, 'arcticdem_06_greenland_northwest'); continue; end

    stripDir = [regionRootDir,'/',regionFolder,'/',region_subdir];
    if exist(stripDir,'dir') ~= 7
        fprintf('%s does not have a subdir "%s", skipping\n', regionFolder, region_subdir);
        continue;
    end

    qcFile = [stripDir,'/','qc.mat'];
    qcLockFile = strrep(qcFile, 'qc.mat', 'qc.mat.lock');

%    if exist(qcFile, 'file')
%        b = load(qcFile);
%        fprintf('%-105s: %5d/%5d (%2.f%%) strip segments have been QC-ed\n', qcFile, nnz(b.flag), length(b.flag), 100*nnz(b.flag)/length(b.flag));
%    end
%    continue

    if exist(qcLockFile, 'file') == 2
        user_id = read_lock(qcLockFile);
        if ~isempty(user_id)
            fprintf(2, 'Region is locked by user "%s", skipping: %s\n', user_id, qcLockFile);
            fprintf(2, 'if this region is locked in error, simply OPEN the lock file in a text editor, CLEAR its contents, and SAVE to unlock\n');
            continue;
        end
    end


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

        fprintf('Appending to %s ... ', qcFile);

        b = {};
        demDname_exist = {};
        if exist(qcFile, 'file') == 2
            b = load(qcFile);
            [demDir_exist,~,~] = cellfun(@fileparts, b.fileNames, 'UniformOutput', false);
            [~,demDname_exist,~] = cellfun(@fileparts, demDir_exist, 'UniformOutput', false);
            demDname_exist = unique(demDname_exist);
        end

        fileNames = {};

        pdir=dir([stripDir,'/*']);
        for j=1:length(pdir);
            %% don't consider this strip if it's already in existing database
            if ~isempty(demDname_exist) && any(strcmp(demDname_exist, pdir(j).name)); continue; end

            demDir = [pdir(j).folder,'/',pdir(j).name];

            finfilecheck=dir([demDir,'/*.fin']);
            if isempty(finfilecheck); continue; end

            fileNamesTemp = dir([demDir,'/*dem_10m_shade_masked.tif']);
            fileNames_add = fullfile({fileNamesTemp.folder},{fileNamesTemp.name});
            %fileNames_add = cellfun(@(x) [x.folder,'/',x.name], {fileNamesTemp}, 'uniformOutput',false);
            fileNames_add = fileNames_add(:);

%            if exist('b','var')
%                [~,IA] = intersect(fileNames_add, b.fileNames);
%                fileNames_add(IA) = [];
%                if isempty(fileNames_add)
%                    continue
%                end
%            end

            fileNames = [fileNames; fileNames_add];
        end

        fprintf('%d new strips\n', length(fileNames));

        if report_number_of_strips_to_append_but_dont_actually_append
            continue;
        end

        if length(fileNames) ~= 0

            flag = zeros(size(fileNames));
            x = cell(size(fileNames));
            y = cell(size(fileNames));

            if ~isempty(b)
                fileNames = [b.fileNames; fileNames];
                flag = [b.flag; flag];
                x = [b.x; x];
                y = [b.y; y];
            end

    %        qcFile = strrep(qcFile, 'qc.mat', 'qc_appended.mat');
            save(qcFile, 'fileNames','x','y','flag');

        end

    else

        fprintf('Creating %s ... ', qcFile);

%        fileNamesTemp = dir([stripDir,'/*/*dem_10m_shade_masked.tif']);
%        fileNames = fullfile({fileNamesTemp.folder},{fileNamesTemp.name});
%        %fileNames = cellfun(@(x) [x.folder,'/',x.name], {fileNamesTemp}, 'uniformOutput',false);
%        fileNames = fileNames(:);

        fileNames = {};

        pdir=dir([stripDir,'/*']);
        for j=1:length(pdir);

            demDir = [pdir(j).folder,'/',pdir(j).name];

            finfilecheck=dir([demDir,'/*.fin']);
            if isempty(finfilecheck); continue; end

            fileNamesTemp = dir([demDir,'/*dem_10m_shade_masked.tif']);
            fileNames_add = fullfile({fileNamesTemp.folder},{fileNamesTemp.name});
            %fileNames_add = cellfun(@(x) [x.folder,'/',x.name], {fileNamesTemp}, 'uniformOutput',false);
            fileNames_add = fileNames_add(:);

            fileNames = [fileNames; fileNames_add];
        end

        fprintf('%d strips\n', length(fileNames))

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


function [user_id] = read_lock(qcLockFile)

    lock_txt = fileread(qcLockFile);
    user_id = strrep(lock_txt, char(10), '');
    user_id = strrep(user_id, '\n', '');

end
