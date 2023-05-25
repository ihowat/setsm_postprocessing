function [z,x,y,z_at_zr_res,I_ref,C] = loadSlopefiltWaterfillArrays(matfile_or_z,x,y,refDemTif,coverTif,varargin)

n = find(strcmpi('skipRefDemLoad',varargin));
if ~isempty(n)
    skipRefDemLoad = true;
else
    skipRefDemLoad = false;
end

m = [];
if isa(matfile_or_z, 'char') || isa(matfile_or_z, 'string')
    if exist(matfile_or_z, 'file')
        m = matfile(matfile_or_z);
    else
        error('matfile_or_z file does not exist: %s', matfile_or_z);
    end
elseif isa(matfile_or_z, 'matlab.io.MatFile')
    m = matfile_or_z;
else
    z = matfile_or_z;
end

if ~isempty(m)
    z = m.z;
    x = m.x;
    y = m.y;
end
dx = x(2)-x(1);

I_ref = readGeotiff(refDemTif,'mapinfoonly');
zr_dx = I_ref.x(2)-I_ref.x(1);
if zr_dx ~= dx
    m_x0 = floor(x(1)  /zr_dx) * zr_dx;
    m_x1 = ceil( x(end)/zr_dx) * zr_dx;
    m_y0 = floor(y(end)/zr_dx) * zr_dx;
    m_y1 = ceil( y(1)  /zr_dx) * zr_dx;
    m_x = m_x0:zr_dx:m_x1;
    m_y = m_y1:-zr_dx:m_y0;
    z_at_zr_res = interp2(x,y(:),z,m_x,m_y(:),'*bilinear');
else
    m_x0 = x(1);
    m_x1 = x(end);
    m_y0 = y(end);
    m_y1 = y(1);
    m_x = x;
    m_y = y;
    z_at_zr_res = z;
end
if ~skipRefDemLoad
    I_ref = getDataFromTileAndNeighbors(refDemTif,m_x0,m_x1,m_y0,m_y1,zr_dx,nan);
end

if isempty(coverTif)
    C = [];
else
    C = readGeotiff(coverTif);
    C.z(C.z == 0) = 80;  % Set NoData (which should be ocean) to water
    C.z = interp2(C.x,C.y(:),C.z,I_ref.x,I_ref.y(:),'*nearest');
end
