function N =refStrip2Mosaic(metaFile,m,dy0,N,varargin)
% refStrip2Mosaic add a reference strip to the mosaic

% read this strip data
[x,y,z,mt,or,c,r] = readStripInTile(metaFile,m.x,m.y);

c=find(c);
r=find(r);

if length(c)*length(r) < 1000; 
    fprintf('%d pixel overlap is too small, skipping\n',sum(c)*sum(r)); 
    m.dtrans=[m.dtrans,[NaN;NaN;NaN]];
    m.rmse=[m.rmse,NaN];
    return; 
end

% check for new coverage
Nsub=N(r,c); %subset N overlap count array

if sum(Nsub(:)==0 & ~isnan(z(:))) < 1000;
    disp('redundant coverage, skipping');
    m.dtrans=[m.dtrans,[NaN;NaN;NaN]];
    m.rmse=[m.rmse,NaN];
    return;
end

% apply registration if given
if length(varargin) == 1
    trans=varargin{1};
    
    x = x + trans(2);
    y = y + trans(3);
    z = z + trans(1);

end

% interpolate the first dem to the same grid
[z,mt,or] = interpolate2grid(x,y,z,mt,or,m.x(1,c),m.y(r,1));

% day number
dy = int16(~isnan(z) .* dy0);

% check if data exists
if any(any(Nsub))
    
    % overlap use feathering
    dx = m.x(1,2)-m.x(1,1);
    buff=10*dx+1;
    W = edgeFeather(Nsub~=0,~isnan(z),buff);

    % make weighted elevation grid
    zsub=m.z(r,c);
    A=  zsub.*W + z.*(1-W);
    n = isnan(zsub) & ~isnan(z);  A(n)= z(n);
    n = ~isnan(zsub) & isnan(z);  A(n)= zsub(n);
    z = A;
    clear A
    
    Nsub = Nsub + uint8(W==0);
    N(r,c)=Nsub;
    clear Nsub
    
    % for the matchtag, just straight combination
    mt = m.mt(r,c) | mt;
    
    % make weighted ortho grid
    or = single(or);
    or(or ==0) = NaN;
    orsub=single(m.or(r,c));
    orsub(orsub==0) = NaN;
    A=  orsub.*W + or.*(1-W);
    
    A( isnan(orsub) & ~isnan(or))=   or(  isnan(orsub) & ~isnan(or));
    A(~isnan(orsub) &  isnan(or))= orsub(~isnan(orsub) &  isnan(or));
    A(isnan(A)) = 0; % convert back to uint16
    or = uint16(A);
    clear A
    
    % make weighted dy grid
    dy = single(dy);
    dy(dy==0) = NaN;
    dysub = single(m.dy(r,c));
    dysub(dysub==0) = NaN;
    
    A=  dysub.*W + dy.*(1-W);
    
    A( isnan(dysub) & ~isnan(dy))=   dy ( isnan(dysub) & ~isnan(dy));
    A(~isnan(dysub) &  isnan(dy))= dysub(~isnan(dysub) &  isnan(dy));
    
    A(isnan(A)) = 0; % convert back to uint16
    dy = uint16(A);
    clear A
    
    clear W
    
else
    N(r,c) = uint8(~isnan(z));
end


m.z(r,c) =  z;
m.mt(r,c) = mt;
m.or(r,c) = or;
m.dy(r,c) = dy;

m.dtrans= [m.dtrans,zeros(3,1)];
m.rmse  = [m.rmse,0];
