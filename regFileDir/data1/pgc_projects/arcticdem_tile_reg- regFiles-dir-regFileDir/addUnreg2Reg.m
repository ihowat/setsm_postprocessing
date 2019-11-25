regFileDir='/data1/pgc_projects/arcticdem_tile_reg';
regFiles=dir([regFileDir,'/*reg_dem.mat']);
regFiles=cellfun( @(x) [regFileDir,'/',x], {regFiles.name},'uniformoutput',0);

i=1;
for i=1:length(regFiles)

regFile=regFiles{i};
unregFile=strrep(regFile,'_reg_','_');

mreg=matfile(regFile,'writable',true);

% check for NaN fields in the registration offsets to indicate missing reg
dtrans = mreg.dtrans
if ~any(sum(~isnan(dtrans)) == 0)
        %none found so skip
        fprintf('%s no NaN dtrans, skipping\n')
        continue
end

% NaN dtrans cols detected, so open files

% load unregfile
munreg = matfile(unregFile);

% load z fields from both files
zreg=mreg.z;
zunreg=munreg.z;

% find pix where reg is missing but unreg has data
N = isnan(zreg) & ~isnan(zunreg);

% check to make sure there are such pix
if any(N(:));

        % add missing unreg pix to reg
        zreg(N) = zunreg(N);
        mreg.z = zreg;

        % load coreg cluster index field from unreg data
        C = munreg.C;
        
        % set registered pix to 1 in coreg cluster, leaving the unreg pix
        % as their cluster #
        C(C ~= 0 & ~N) = 1;

        % put field into reg file
        mreg.C = C;
end

end
