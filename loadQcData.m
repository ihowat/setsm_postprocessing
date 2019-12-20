function meta = loadQcData(dbasefile)
%loadQcData - Loads QC.mat file data from folders with DEMs in the meta db
%   into the meta db structure


% load database file into mat file object
meta=load(dbasefile);

regionDir = regexprep(meta.f, '[\d\w]*.txt$', '');
%disp(regionDir);
uniqueRegionDir=unique(regionDir);
%disp(uniqueRegionDir);
i=1;
qc.fileNames=[];
qc.flag=[];
qc.x =[];
qc.y=[];
for i=1:length(uniqueRegionDir)

    if exist([uniqueRegionDir{i},'/qc.mat'],'file')

        qci=load([uniqueRegionDir{i},'/qc.mat']);

        qc.fileNames=[qc.fileNames(:);qci.fileNames(:)];
        qc.flag=[qc.flag(:);qci.flag(:)];
        qc.x=[qc.x(:);qci.x(:)];
        qc.y=[qc.y(:);qci.y(:)];

        clear qci;
    end

end

if ~isempty(qc.fileNames)
    
    qc.fileNames=strrep(qc.fileNames,'V:','/mnt');
    qc.fileNames=strrep(qc.fileNames,'\','/');
    
    [~,IA,IB]=intersect( strrep(meta.f,'_meta.txt',''), strrep(qc.fileNames,'_dem_10m_shade_masked.tif',''));

    fprintf('%d of %d strips with qc data\n',length(IA),length(meta.f))

    meta.qc = zeros(size(meta.f),'uint8');
    meta.maskPolyx = cell(size(meta.f));
    meta.maskPolyy = cell(size(meta.f));

    meta.qc(IA) = qc.flag(IB);
    meta.maskPolyx(IA) = qc.x(IB);
    meta.maskPolyy(IA) = qc.y(IB);

    % apply filter
    n=meta.qc >=1 & meta.qc <=3;
    %n = meta.qc >=1 ;

    meta = structfun(@(x) ( x(n,:) ), meta, 'UniformOutput', false);

    % qc 3 to 2 so that they are treated equally in the mosaicker
    meta.qc(meta.qc == 3) = 2;

end


