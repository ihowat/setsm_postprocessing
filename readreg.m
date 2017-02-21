function [N,mn,md,trans,pctile]=readreg(f)



fid=fopen(f);
c=textscan(fid,'%s','delimiter','\n');
fclose(fid);
c=c{1};

astr='# GCPs=';
r=~cellfun(@isempty,strfind(c,astr));
N = str2double(strrep(c{r},astr,''));


astr='All GCP DZ AVG=';
r=~cellfun(@isempty,strfind(c,astr));
mn = str2double(strrep(c{r},astr,''));


astr='All GCP DZ MED=';
r=~cellfun(@isempty,strfind(c,astr));
md = str2double(strrep(c{r},astr,''));


astr='Translation vector (dz,dx,dy)(meters)=';
r=~cellfun(@isempty,strfind(c,astr));
r= strrep(c{r},astr,'');
trans = sscanf(r,'%f, %f, %f');
trans = trans(:)';

% p=50:5:100;
% j=1;
% for j=1:length(p)
%     astr=[num2str(p(j)),'th percentile='];
%     r=find(~cellfun(@isempty,strfind(c,astr)));
%     pctile(1,j)= str2double(strrep(c{r},astr,''));
% end

astr='All=';
r=~cellfun(@isempty,strfind(c,astr));
pctile = str2num(strrep(c{r},astr,''));