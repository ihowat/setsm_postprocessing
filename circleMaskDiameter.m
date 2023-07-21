function mask = circleMaskDiameter(diameter)

if mod(diameter, 2) == 0
   diameter = diameter - 1;
end
radius = (diameter - 1) / 2;

dim = radius * 2 + 1;

xdim = dim;
ydim = dim;

xc = radius + 1;
yc = radius + 1;

[xx,yy] = meshgrid(1:xdim, 1:ydim);
mask = false(ydim, xdim);
mask = mask | hypot(xx - xc, yy - yc) <= radius;
