function [N,mn,md,trans,pctile]=readisreg(f)



fid=fopen(f);
c=textscan(fid,'%s','delimiter','\n');
fclose(fid);
c=c{1};

astr='# GCPs=';
r=~cellfun(@isempty,strfind(c,astr));
N = str2double(strrep(c{r},astr,''));


astr='Mean Vertical Residual (m)=';
r=~cellfun(@isempty,strfind(c,astr));
mn = str2double(strrep(c{r},astr,''));


astr='Median Vertical Residual (m)=';
r=~cellfun(@isempty,strfind(c,astr));
md = str2double(strrep(c{r},astr,''));

astr='Translation Vector (dz,dx,dy)(m)=';
r=~cellfun(@isempty,strfind(c,astr));
trans = str2num(strrep(c{r},astr,''));

p=50:5:100;
j=1;
for j=1:length(p)
    astr=[num2str(p(j)),'th percentile='];
    r=find(~cellfun(@isempty,strfind(c,astr)));
    pctile(1,j)= str2double(strrep(c{r},astr,''));
end
