function [z,mt,or,dy] = refStrip2Mosaic(metaFile,x,y,z,mt,or,dy,dy0)
% refStrip2Mosaic add a reference strip to the mosaic

[d,m,o,c,r] = readStrip(metaFile,x,y);

if sum(c)*sum(r) < 100; disp('too little overlap'); return; end

% interpolate the first dem to the same grid
zi = interp2(d.x,d.y',double(d.z),...
    x(c),y(r),'*linear');

[mi,oi] = interpolate2grid(m,o,x(c),y(r)',[0;0;0]);

zsub=z(r,c);
mtsub=mt(r,c);
orsub=or(r,c);
dysub=dy(r,c);

dyi = ~isnan(zi) .* dy0;

zi(isnan(zi) & ~isnan(zsub)) = zsub(isnan(zi) & ~isnan(zsub));

oi(oi == 0 & orsub ~= 0) = orsub(oi == 0 & orsub ~= 0);

dyi(dyi == 0 & dysub ~= 0) = dysub(dyi == 0 & dysub ~= 0);

z(r,c) = zi;
mt(r,c) = mi | mtsub;
or(r,c) = oi;
dy(r,c) = int16(dyi);
