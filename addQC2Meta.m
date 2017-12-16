function m = addQC2Meta(m,changePath)

m.baseName = strrep(m.f,'_8m_meta.txt','');

uniqueRegionDir=unique(m.region);
i=1;
qc.fileNames=[];
qc.flag=[];
qc.x =[];
qc.y=[];
for i=1:length(uniqueRegionDir)
    
    if exist([uniqueRegionDir{i},'/qc.mat'],'file')
        
        qci=load([uniqueRegionDir{i},'/qc.mat']);
                
        n = ~cellfun(@iscell, qci.x) | ~cellfun(@iscell, qci.y);
        if any(n)
            warning('non cell mask entries detected, correcting')
            n=find(n);
            j=1;
            for j=1:length(n)
                qci.x(n(j)) = {{qci.x{n(j)}'}};
                qci.y(n(j)) = {{qci.y{n(j)}'}};
            end
            
            save([uniqueRegionDir{i},'/qc.mat'],'-struct','qci')
            
        end
        
        qc.fileNames=[qc.fileNames(:);qci.fileNames(:)];
        qc.flag=[qc.flag(:);qci.flag(:)];
        qc.x=[qc.x(:);qci.x(:)];
        qc.y=[qc.y(:);qci.y(:)];
        
        clear qci;
    end
    
end

if ~isempty(changePath)
    qc.fileNames= strrep(qc.fileNames,'/data4',changePath);
end

if ~isempty(qc.fileNames)
    
    [~,IA,IB]=intersect( strrep(m.f,'_meta.txt',''), strrep(qc.fileNames,'_dem_browse.tif',''));
    
    if isempty(IA)
        [~,IA,IB]=intersect( strrep(m.f,'_meta.txt',''), strrep(qc.fileNames,'_dem.tif',''));
    end
   
    m.qc = zeros(size(m.f),'uint8');
    m.maskPolyx = cell(size(m.f));
    m.maskPolyy = cell(size(m.f));
    
    m.qc(IA) = qc.flag(IB);
    m.maskPolyx(IA) = qc.x(IB);
    m.maskPolyy(IA) = qc.y(IB);
    
end