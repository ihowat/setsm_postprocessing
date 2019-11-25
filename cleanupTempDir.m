function cleanupTempDir()

global TEMPDIR;
if ~isempty(TEMPDIR)
    if exist(TEMPDIR, 'dir') == 7
        fprintf('Removing temporary directory for this run: %s\n', TEMPDIR);
        rmdir(TEMPDIR);
    else
        fprintf('Temporary directory was removed prematurely? %s\n', TEMPDIR);
    end
end

TEMPDIR = [];
