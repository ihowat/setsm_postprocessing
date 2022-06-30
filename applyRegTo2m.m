function applyRegTo2m(fileNames,dir2m)
% applyDzfitTo2m apply the dzfit field from 10m to its 2m quad tiles
%
% applyDzfitTo2m(fileNames,dir2m) with the 10m files containing the dzfit
% field and the path in dir2m containing the corresponding 2m quad tiles.
% Insert '.' if in same path. The 10m files must have a name format 
% *10m_reg.mat with the corresponding 2m having a format *2m_r_c_reg.mat


if ~iscell(fileNames)
    fileNames={fileNames};
end

dir2m_arg = dir2m;


% loop through files (can be done in parallel)
i=1;
for i=1:length(fileNames)
    
    % get 10m file name
    fileName10m=fileNames{i};

    fprintf('Working on 10m tile: %s\n', fileName10m)
    
    [tiledir,filename10m_base,filename10m_ext]= fileparts(fileName10m);
    fileName10mName = [filename10m_base,filename10m_ext];

    if isempty(dir2m_arg)
        dir2m = tiledir;
    end
    fileNames2m = {[dir2m,'/',strrep(fileName10mName,'10m_reg.mat','1_1_2m.mat')],...
                   [dir2m,'/',strrep(fileName10mName,'10m_reg.mat','1_2_2m.mat')],...
                   [dir2m,'/',strrep(fileName10mName,'10m_reg.mat','2_1_2m.mat')],...
                   [dir2m,'/',strrep(fileName10mName,'10m_reg.mat','2_2_2m.mat')]};
    
    n = cellfun( @exist, fileNames2m);
    
    if ~any(n)
        warning('ERROR: no 2m files found for %s',fileName10m)
        continue
    end
    
    fileNames2m(~n) = [];

    m10 = matfile(fileName10m);

    % get the size of 10m z
    s=whos(m10);
    sz = s(strcmp({s.name},'z')).size;

    if any(strcmp(fields(m10),'dzfit')) && ~isempty(m10.dzfit)
        % resize dzfit to 10m z
        dzfit0 = imresize(m10.dzfit,sz);
    else
        dzfit0 = [];
    end
    
    
    j=1;
    for j=1:length(fileNames2m)

        fileName2m_unreg = fileNames2m{j};
        fileName2m_reg = strrep(fileName2m_unreg,'.mat','_reg.mat');

        fprintf('Working on 2m tile: %s\n', fileName2m_unreg)

        applyRegistration(fileName2m_unreg,[],'demMatFileWithRegToCopy',fileName10m);
        if ~exist(fileName2m_reg)
            error('failed to create 2m reg tile file for %s',fileName2m_unreg)
        end

        fprintf('2m reg file was created: %s\n', fileName2m_reg)

        if ~isempty(dzfit0)
            fprintf('Applying dzfit from 10m to 2m tile\n')

            m2 = matfile(fileName2m_reg);

            % interpolate 10m full tile to 2m quad grid
            dzfit = interp2(m10.x,m10.y,dzfit0,m2.x,m2.y,'*bilinear');

            % apply to z
            m2.Properties.Writable = true;
            m2.z = m2.z - dzfit;
            m2.dzfitApplied = true;

            % downsample dzfit and add to mat file
            m2.dzfit = imresize(dzfit,0.01);
        end
    end
    
end
    
   
    
    
        
    
    
