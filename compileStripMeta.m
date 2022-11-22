function [f,creationDate,stripDate,res,x,y,avg_rmse,med_rmse,max_rmse,Nscenes,reg]=compileStripMeta(fdir)
meta_str='meta';
reg_str='reg';
disableReg=true; 

f=dir([fdir,'/*',meta_str,'.txt']);
f={f.name}';

creationDate=nan(size(f));
x=cell(size(f));
y=cell(size(f));

res=nan(size(f));
avg_rmse=nan(size(f));
med_rmse=nan(size(f));
max_rmse=nan(size(f));
Nscenes=nan(size(f));

stripDate=arrayfun(@(x) datenum(parsePairnameDatestring(convertCharsToStrings(x)),'yyyymmdd'), f, 'uniformoutput',0);

 reg.datasets       = cell(size(f));
 reg.Ngcps          = nan(size(f));
 reg.NgcpsRock      = nan(size(f));
 reg.NgcpsIce       = nan(size(f));
 reg.quads          = nan(size(f));
 reg.avgIceGCPdt    = nan(size(f));
 reg.trans          = nan(size(f,1),3);
 reg.avgdz          = nan(size(f));
 reg.avgdz_rock     = nan(size(f));
 reg.avgdz_ice      = nan(size(f));
 reg.stddz          = nan(size(f));
 reg.stddz_rock     = nan(size(f));
 reg.stddz_ice      = nan(size(f));
 reg.dzpc           = nan(size(f,1),11);
 reg.dzpc_rock      = nan(size(f,1),11);
 reg.dzpc_ice       = nan(size(f,1),11);


for i=1:length(f)

    c=textread([fdir,'/',f{i}],'%s','delimiter','\n');
    
    r=find(~cellfun(@isempty,strfind(c,'Strip creation date:')));
    creationDate(i)=datenum(deblank(strrep(c{r},'Strip creation date:','')));
    
    info=imfinfo(strrep([fdir,'/',f{i}],'_meta.txt','_dem.tif'));
    res(i)=info.ModelPixelScaleTag(1);
    
    r=find(~cellfun(@isempty,strfind(c,'X:')));
    x{i}=str2num(strrep(c{r},'X:',''));
    
    r=find(~cellfun(@isempty,strfind(c,'Y:')));
    y{i}=str2num(strrep(c{r},'Y:',''));
    
    
    r=find(~cellfun(@isempty,strfind(c,'Mosaicking Alignment Statistics')));
    r=r+3;
    A=[];
    while r
        if isempty(c{r}); break; end
        A=[A;sscanf(c{r},'%*s %f %f %f %f')'];
        r=r+1;
    end
    
    clear c
    
    if isempty(A);
        Nscenes(i)=1;
    else
        A(:,1) = round(A(:,1), 2);
        
        avg_rmse(i)=nanmean(A(:,1));
        med_rmse(i)=nanmedian(A(:,1));
        max_rmse(i)=nanmax(A(:,1));
        Nscenes(i)=size(A,1)+1;
        
    end
   
   
    regfile=strrep([fdir,'/',f{i}],meta_str,reg_str);
    if exist(regfile,'file') && ~disableReg
        
        c=textread(regfile,'%s','delimiter','\n');
        
        astr='Registration Dataset 1 Name:';
        r=find(~cellfun(@isempty,strfind(c,astr)));
        if length(r) > 1; disp('multiple datasets'); end
        reg.datasets{i} = strrep(c{r},astr,'');

        astr='# GCPs=';
        r=find(~cellfun(@isempty,strfind(c,astr)));
        reg.Ngcps(i) = str2num(strrep(c{r},astr,''));
        
        astr='# Rock GCPs=';
        r=find(~cellfun(@isempty,strfind(c,astr)));
        if (~exist('r'))
            reg.NgcpsRock(i) = str2num(strrep(c{r},astr,''));
        end
        
        astr='# Ice GCPs=';
        r=find(~cellfun(@isempty,strfind(c,astr)));
        if (~exist('r'))
            reg.NgcpsIce(i) = str2num(strrep(c{r},astr,''));
        end
        
        astr='Quadrants with GCPs=';
        r=find(~cellfun(@isempty,strfind(c,astr)));
        if (~exist('r'))
            reg.quads(i) = str2num(strrep(c{r},astr,''));
        end
        
        if  reg.NgcpsIce(i) > 0;
            astr='Average Ice GCP time offset=';
            r=find(~cellfun(@isempty,strfind(c,astr)));
            reg.avgIceGCPdt(i)= str2num(strrep(c{r},astr,''));
        end
        
        astr='Translation vector (dz,dx,dy)(meters)=';
        r=find(~cellfun(@isempty,strfind(c,astr)));
        reg.trans(i,:)= str2num(strrep(c{r},astr,''));
        
        astr='All GCP DZ AVG=';
        r=find(~cellfun(@isempty,strfind(c,astr)));
        if (~exist('r'))
            reg.avgdz(i)= str2num(strrep(c{r},astr,''));
        end
        
        astr='Rock GCP DZ AVG=';
        r=find(~cellfun(@isempty,strfind(c,astr)));
        if (~exist('r'))
            reg.avgdz_rock(i)= str2num(strrep(c{r},astr,''));
        end
        
        astr='Ice GCP DZ AVG=';
        r=find(~cellfun(@isempty,strfind(c,astr)));
        if (~exist('r'))
            reg.avgdz_ice(i)= str2num(strrep(c{r},astr,''));
        end
        
        astr='All GCP DZ MED=';
        r=find(~cellfun(@isempty,strfind(c,astr)));
        if (~exist('r'))
            reg.avgdz(i)= str2num(strrep(c{r},astr,''));
        end
        
        astr='Rock GCP DZ MED=';
        r=find(~cellfun(@isempty,strfind(c,astr)));
        if (~exist('r'))
            reg.avgdz_rock(i)= str2num(strrep(c{r},astr,''));
        end
        
        astr='Ice GCP DZ MED=';
        r=find(~cellfun(@isempty,strfind(c,astr)));
        if (~exist('r'))
            reg.avgdz_ice(i)= str2num(strrep(c{r},astr,''));
        end
        
        astr='All GCP DZ STD=';
        r=find(~cellfun(@isempty,strfind(c,astr)));
        if (~exist('r'))
            reg.stddz(i)= str2num(strrep(c{r},astr,''));
        end
        
        astr='Rock GCP DZ STD=';
        r=find(~cellfun(@isempty,strfind(c,astr)));
        if (~exist('r'))
            reg.stddz_rock(i)= str2num(strrep(c{r},astr,''));
        end
        
        astr='Ice GCP DZ STD=';
        r=find(~cellfun(@isempty,strfind(c,astr)));
        if (~exist('r'))
            reg.stddz_ice(i)= str2num(strrep(c{r},astr,''));
        end
        
        astr='All=';
        r=find(~cellfun(@isempty,strfind(c,astr)));
        if (~exist('r'))
            reg.dzpc(i,:) = str2num(strrep(c{r},astr,''));
        end
        
        astr='Rock=';
        r=find(~cellfun(@isempty,strfind(c,astr)));
        if (~exist('r'))
            reg.dzpc_rock(i,:)= str2num(strrep(c{r},astr,''));
        end
        
        astr='Ice=';
        r=find(~cellfun(@isempty,strfind(c,astr)));
        if (~exist('r'))
            reg.dzpc_ice(i,:) = str2num(strrep(c{r},astr,''));
        end
        
    end
    
end


 
