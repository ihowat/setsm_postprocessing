function fit2is2(fileNames,is2dir)
%fit2is2 fit DEM tiles to icesat2 points using a quadratic surface
%
%fit2is2(fileNames,is2dir) is a cell array of tile file names and is2dir is
% directory of the icesat 2 matlab files. The tile files must be of the
% form rr_cc_10m_reg.mat and the icesat2 files must be rr_cc_is2.mat . Adds
% an offset array dzfit to the tile file that is 1/100 the size of z. This
% is imresized to the dims of z and applied as z - dzfit. 

% make sure input is cell array - in case just one file name is passed as a
% string
if ~iscell(fileNames)
    fileNames={fileNames};
end

% loop through files (can be done in parallel)
i=1;
for i=1:length(fileNames)
    
    fileName = fileNames{i};
    
    fprintf('fitting %s\n',fileName)
    
    is2name=[is2dir,'/',strrep(fileName,'_10m_reg.mat','_is2.mat')];
    
    if ~exist(is2name,'file')
        warning('%s not found, skipping',is2name)
    end
    
    m=matfile(fileName);
    
    is2 =load(is2name);
    
    % get the size of z
    s=whos(m);
    if any(strcmp({s.name},'land'))
        [dzfit,sf]=fitDEM2gcps(m.x,m.y,m.z,is2.x,is2.y,is2.z,'resizeFactor',0.01,'landMask',m.land);
    else
        [dzfit,sf]=fitDEM2gcps(m.x,m.y,m.z,is2.x,is2.y,is2.z,'resizeFactor',0.01);
    end
    
    if isempty(dzfit)
        continue
    end
    
    % add downscaled surface to the file
    m.Properties.Writable = true;
    m.dzfit = dzfit;
    m.dzfit = sf;
    
    % get the size of z
    sz = s(strcmp({s.name},'z')).size;
    
    % resize dzfit to z
    dzfit = imresize(dzfit,sz);
    
    % apply offset
    m.z = m.z - dzfit;
    
    clear dzfit is2 fileName sf
end








