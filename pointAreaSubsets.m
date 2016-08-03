function [xsub,ysub,zsub]=pointAreaSubsets(xp,yp,x,y,m,dd,N)

% initialize coordinate subset cells
xsub=cell(length(xp),1);
ysub=cell(length(xp),1);
zsub=cell(length(xp),1);

% subsetting loop
j=1;
for j=1:length(xp)
    minCol=find(x >=  xp(j) - dd,1,'first');
    maxCol=find(x <=  xp(j) + dd,1,'last');
    
    minRow=find(y <=  yp(j) + dd,1,'first');
    maxRow=find(y >=  yp(j) - dd,1,'last');
    
    xsub{j}=x(minCol:maxCol);
    ysub{j}=y(minRow:maxRow);
    zsub{j}=m.z(minRow:maxRow,minCol:maxCol);
    zsub{j}(~N(minRow:maxRow,minCol:maxCol)) = NaN;
end

