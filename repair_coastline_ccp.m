load('PGC_Imagery_Mosaic_Tiles_Arctic.mat');
tns = {'58_23','58_24','59_23','59_24','60_22','60_23','60_24'};
for j=1:length(tns)
    i=strcmp(I,tns{j});
    coastline{i}=[x0(i),x0(i),x1(i),x1(i),x0(i);...
        y0(i),y1(i),y1(i),y0(i),y0(i)];
end
save('PGC_Imagery_Mosaic_Tiles_Arctic_coastrepair.mat','I','x0','x1','y0','y1','coastline','-v7.3');