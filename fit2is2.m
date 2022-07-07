function fit2is2(tileFileName,is2FileName,varargin)
%fit2is2 fit DEM tiles to icesat2 points using a quadratic surface
%
%fit2is2(fileNames,is2dir) is a cell array of tile file names and is2dir is
% directory of the icesat 2 matlab files. The tile files must be of the
% form rr_cc_10m_reg.mat and the icesat2 files must be rr_cc_is2.mat . Adds
% an offset array dzfit to the tile file that is 1/100 the size of z. This
% is imresized to the dims of z and applied as z - dzfit.
m=matfile(tileFileName);
s=whos(m);

% get the size of z
sz = s(strcmp({s.name},'z')).size;

dzfitMinPoints = [];
n = find(strcmpi('dzfitMinPoints', varargin));
if ~isempty(n)
    dzfitMinPoints = varargin{n+1};
end

if any(strcmp(fields(m),'dzfitApplied'))
    if m.dzfitApplied==1
        if ~any(strcmp(varargin,'overwrite'))
            fprintf('dzfit already applied, returning\n')
            return
        elseif any(strcmpi(fields(m),'dzfit'))
            fprintf('undoing dzfit\n')
            m.Properties.Writable = true;
            dzfit = imresize(m.dzfit,sz);
            m.z = m.z + dzfit;
            m.dzfitApplied=false;
            fprintf('proceeding to recalc and reapply dzfit\n')
        end
    end
end
%load is2 data
is2 =load(is2FileName,'x','y','z');

try
    %load DEM
    load(tileFileName,'x','y','z','land');
catch
    warning('cant read %s\n',tileFileName)
    return
end

% remove is2 points that fall outside of the (2m) tile
n = is2.x >= min(x) & is2.x <= max(x) & is2.y >= min(y) & is2.y <= max(y);
is2 = structfun( @(x) x(n), is2, 'uniformoutput', 0);

% call the surfacef fitter with or without land mask if it exists
if any(strcmp({s.name},'land'))
    [dzfit,sf]=fitDEM2gcps(x,y,z,is2.x,is2.y,is2.z,'resizeFactor',0.01,'landMask',land,'minPoints',dzfitMinPoints);
else
    [dzfit,sf]=fitDEM2gcps(x,y,z,is2.x,is2.y,is2.z,'resizeFactor',0.01,'minPoints',dzfitMinPoints);
end

if isempty(dzfit)
    fprintf('dzfit is empty and will not be applied\n')
    return
else
    fprintf('applying dzfit\n')
end

% add downscaled surface to the file
m.Properties.Writable = true;
m.dzfit = dzfit;
m.sf = sf;

% resize dzfit to z
dzfit = imresize(dzfit,sz);

% apply offset
m.z = m.z - dzfit;
m.dzfitApplied = true;








