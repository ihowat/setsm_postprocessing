function tempdir = getTempDir(local_dir)

global TEMPDIR;
if ~isempty(TEMPDIR)
    tempdir = TEMPDIR;
    return;
end

global TEMPDIR_CLEANUP;
TEMPDIR_CLEANUP = onCleanup(@() cleanupTempDir());

if ~exist('local_dir', 'var')
    if ispc
        local_dir = getenv('USERPROFILE');
    elseif isunix
        local_dir = '/local/';
        if exist(local_dir, 'dir') ~= 7
            local_dir = fullfile(getenv('HOME'), 'scratch');
        end
    end
    local_dir = fullfile(local_dir, 'setsm_postprocessing_temp');
    
    if exist(local_dir, 'dir') ~= 7
        fprintf('Creating new root directory for temporary files: %s\n', local_dir);
        mkdir(local_dir);
    end
    
elseif exist(local_dir, 'dir') ~= 7
    error('Specified `local_dir` directory does not exist: %s', local_dir);
end

tempdir = tempname(local_dir);

if exist(tempdir, 'dir') ~= 7
    fprintf('Creating temporary directory for this run: %s\n', tempdir);
    mkdir(tempdir);
end

TEMPDIR = tempdir;
