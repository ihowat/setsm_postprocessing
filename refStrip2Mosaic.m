function N =refStrip2Mosaic(metaFile,m,dy0,N,varargin)
% refStrip2Mosaic add a reference strip to the mosaic



% read this strip data
[x,y,z,mt,or,c,r] = readStripInTile(metaFile,m.x,m.y);

c=find(c);
r=find(r);

if length(c)*length(r) < 100; 
    fprintf('%d pixel overlap is too small, skipping\n',sum(c)*sum(r)); 
    m.dtrans=[m.dtrans,[NaN;NaN;NaN]];
    m.rmse=[m.rmse,NaN];
    return; 
end

% apply registration if given
if length(varargin) == 1
    trans=varargin{1};
    
    x = x - trans(2);
    y = y - trans(3);
    z = z - trans(1);
    
    
end

% interpolate the first dem to the same grid
[z,mt,or] = interpolate2grid(x,y,z,mt,or,m.x(1,c),m.y(r,1));

dy = int16(~isnan(z) .* dy0);

if any(any(N(r,c)))
    
    Nsub=N(r,c);
    Nsub = Nsub + uint8(~isnan(z));
    N(r,c) = Nsub;
    clear Nsub
    
    zsub=m.z(r,c);
    n= isnan(z) &  ~isnan(zsub);
    z(n)=zsub(n);
    clear zsub

    orsub=m.or(r,c);
    n= or==0 & orsub ~= 0;
    or(n)=orsub(n);
    clear orsub
    
    dysub=m.dy(r,c);
    n= dy==0 & dysub ~= 0;
    dy(n)=dysub(n);
    clear dysub
    
    mtsub=m.mt(r,c);
    mt = mt | mtsub;
    clear mtsub
    
else
    N(r,c) = uint8(~isnan(z));
end


m.z(r,c) =  z;
m.mt(r,c) = mt;
m.or(r,c) = or;
m.dy(r,c) = dy;

m.dtrans= [m.dtrans,zeros(3,1)];
m.rmse  = [m.rmse,0];    