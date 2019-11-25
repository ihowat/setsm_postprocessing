function QCStrips(stripDir,varargin)


outFile=[stripDir,'/qc.mat'];

%system(['rm ',stripDir,'/._*']);

fileNames=dir([stripDir,'/*dem_browse.tif']);

fileNames=cellfun( @(x) [stripDir,'/',x], {fileNames.name},'uniformOutput',false);
fileNames=fileNames(:);

flag=zeros(size(fileNames));
x=cell(size(fileNames));
y=cell(size(fileNames));

if exist(outFile,'file')
 
    a=load(outFile);
    
    [~,IA,IB] = intersect(fileNames,a.fileNames);
    
    flag(IA) = a.flag(IB);
    x(IA) = a.x(IB);
    y(IA) = a.y(IB);
    
    clear a
    
end


if ~isempty(varargin)
    
    switch lower(varargin{1})
        
        case 'redolast'
            
            n=find(flag ~=0,1,'last');
            
            flag(n)=0;
            x{n}=[];
            y{n}=[];
    end
end


for i=1:length(fileNames)
    
    if flag(i) ~= 0; continue; end
    
    fprintf('%d of %d\n',sum(flag ~= 0) + 1,length(fileNames))
    
    [x{i},y{i},flag(i)]=imedit(fileNames{i});
    
    save(outFile,'fileNames','x','y','flag');
    %eval(['! chmod 664 ',outFile]);
    
end



function [x,y,flag]=imedit(fileName)

I=readGeotiff(fileName);
%Ior=readGeotiff(strrep(fileName,'_dem_browse.tif','_ortho_browse.tif'));

% %%
% landclass = subsetGimpIceMask(...
%     min(I.x)-30,max(I.x)+30,min(I.y)-30,max(I.y)+30);
% 
% if ~isempty(landclass.ocean)
%     
%     if any(landclass.ocean(:))
%         
%         
%         I.ocean = interp2(landclass.x,landclass.y(:),single(landclass.ocean),...
%             I.x,I.y(:),'*nearest');
%         
%         I.z(logical(I.ocean))=0;
%     end
% end
%%
imagesc(I.x,I.y,I.z,'alphadata',I.z ~= 0)
set(gca,'color','r')
axis xy  equal tight;
colormap gray;
hold on;

% set(gcf,'units','normalized');
% set(gcf,'position',[0.01,0.01,.35,.9])


flag=[];
x=cell(1);
y=cell(1);

i=1;
while i
    fprintf('%s\n',fileName)
    flag=input('Enter quality flag: 1=good, 2=partial, 3=manual edit, 4=poor, 5=bad/error\n');
    
    if ~isempty(flag)
        if isnumeric(flag)
            if flag == 1 || flag == 2 || flag == 3 || flag == 4 || flag == 5
                break;
            end
        end
    end
    
    if iscell(flag); flag=flag{1}; end
    
    fprintf('%d not recogized, try again\n',flag);

end

if flag == 1 || flag == 2; clf; return; end

if flag == 4 || flag == 5
    x{1} = [min(I.x);min(I.x);max(I.x);max(I.x);min(I.x)];
    y{1} = [min(I.y);max(I.y);max(I.y);min(I.y);min(I.y)];
    clf
    return
end

fprintf('entering manual edit mode\n')

while i
    
    [~,x{i},y{i}] =  roipoly;
    
    plot(x{i},y{i},'r')
    
    
    while i
        s=input('continue editing this image? (y/n)\n','s');
    
        if ~strcmpi(s,'y') && ~strcmpi(s,'n')
            fprintf('%s not recogized, try again\n',s);
        else
            break
        end
    end
    
    if strcmpi(s,'n'); break; end
    
    i=i+1;
end

clf

        
    