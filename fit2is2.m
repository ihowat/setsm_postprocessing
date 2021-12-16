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

if any(strcmp(fields(m),'dzfitApplied'))
    if m.dzfitApplied==1 && ~any(strcmp(varargin,'overwrite'))
        fprintf('dzfit already applied, retutning\n')
        return
    end
end
%load is2 data
is2 =load(is2FileName,'x','y','z');

% call the surfacef fitter with or without land mask if it exists
if any(strcmp({s.name},'land'))
    [dzfit,sf]=fitDEM2gcps(m.x,m.y,m.z,is2.x,is2.y,is2.z,'resizeFactor',0.01,'landMask',m.land);
else
    [dzfit,sf]=fitDEM2gcps(m.x,m.y,m.z,is2.x,is2.y,is2.z,'resizeFactor',0.01);
end

if isempty(dzfit)
    return
end

% add downscaled surface to the file
m.Properties.Writable = true;
m.dzfit = dzfit;
m.sf = sf;

% get the size of z
sz = s(strcmp({s.name},'z')).size;

% resize dzfit to z
dzfit = imresize(dzfit,sz);

% apply offset
m.z = m.z - dzfit;
m.dzfitApplied = true;








