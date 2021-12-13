function applyDzTo2m(fileNames,dir2m)
% applyDzfitTo2m apply the dzfit field from 10m to its 2m quad tiles
%
% applyDzfitTo2m(fileNames,dir2m) with the 10m files containing the dzfit
% field and the path in dir2m containing the corresponding 2m quad tiles.
% Insert '.' if in same path. The 10m files must have a name format 
% *10m_reg.mat with the corresponding 2m having a format *2m_r_c_reg.mat


if ~iscell(fileNames)
    fileNames={fileNames};
end


% loop through files (can be done in parallel)
i=1;
for i=1:length(fileNames)
    
    % get 10m file name
    fileName=fileNames{i};
    
    [~,name]= fileparts(fileName);
    
    fileNames2m = {[dir2m,'/',strrep(name,'10m_reg','1_1_2m_reg.mat')],...
                   [dir2m,'/',strrep(name,'10m_reg','1_2_2m_reg.mat')],...
                   [dir2m,'/',strrep(name,'10m_reg','2_1_2m_reg.mat')],...
                   [dir2m,'/',strrep(name,'10m_reg','2_2_2m_reg.mat')]};
    
    n = cellfun( @exist, fileNames2m);
    
    if ~any(n)
        warning('no 2m files found for %s',fileName)
        continue
    end
    
    fileNames2m(~n) = [];

    m10 = matfile(fileName10m);
    
    % get the size of 10m z
    s=whos(m10);
    sz = s(strcmp({s.name},'z')).size;
    
    % resize dzfit to 10m z
    dzfit0 = imresize(m10.dzfit,sz);
    
    % resize dz0 to 10m z
    dz0 = imresize(m10.dz0,sz);
    
    
    j=1;
    for j=1:length(fileNames2m)

        m2 = matfile(fileNames2m{j});
        
        % interpolate 10m full tile to 2m quad grid
        dzfit = interp2(m10.x,m10.y,dzfit0,m2.x,m2.y,'*bilinear');
        dz0 = interp2(m10.x,m10.y,dzfit0,m2.x,m2.y,'*bilinear');
        
        % apply to z
        m2.Properties.Writable = true;
        m2.z = m2.z - dzfit;
        m2.z = m2.z - dz0;
        
        % downsample dzfit and add to mat file
        m2.dzfit = imresize(dzfit,0.01);
        m2.dz0 = imresize(dz0,0.01);
    end
    
end
    
   
    
    
        
    
    
    