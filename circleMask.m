function mask = circleMask(radius)

dim = radius * 2 + 1;
xdim = dim;
ydim = dim;

xc = radius + 1;
yc = radius + 1;

[xx,yy] = meshgrid(1:xdim, 1:ydim);
mask = false(ydim, xdim);
mask = mask | hypot(xx - xc, yy - yc) <= radius;
