function mask = circleFeatherMaskRadius(radius, featherStartRadius)

dim = radius * 2 + 1;
xdim = dim;
ydim = dim;

xc = radius + 1;
yc = radius + 1;

[xx,yy] = meshgrid(1:xdim, 1:ydim);
mask = 1 - min(1, max(0, (hypot(xx - xc, yy - yc) - featherStartRadius) / (radius - featherStartRadius)));
