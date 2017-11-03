function [meta,m,N,C] = regStrips2Tile(meta,m,N,C,mergeMethodReg,dy0,minNewPixels)
% regStrips2Tile add registered strips to tile
%
% [meta,m,N,C] = regStrips2Tile(meta,m,N,C,,mergeMethodReg) 
% determines if registration data exist for strips in meta.f, either in the 
% metastruct (fields dtrans and rmse) or in *_reg.txt files, and passes the
% the strips to the mosaicker. dy0 is the zero day for the image date
% (ArcticDEM uses Y2K)
%
% Ian Howat, ihowat@gmail.com, Ohio State

%% if dtrans is not in the meta, look for reg files.
if ~isfield(meta,'dtrans')
    
    regfiles=strrep(meta.f,'meta.txt','reg.txt');
    nreg= find( cellfun( @exist, regfiles) == 2);
    
    meta.dtrans=nan(size(meta.f,1),3);
    meta.rmse=nan(size(meta.f));
    
    if ~isempty(nreg)
        regfiles=regfiles(nreg);
        
        [~,~,~,dtrans,pctile]=cellfun( @readreg,regfiles,'uniformoutput',0);
        pctile=cell2mat(pctile);
        dtrans=cell2mat(dtrans);
        p75=pctile(:,6);
        
        meta.dtrans(nreg,:)=dtrans;
        meta.rmse(nreg)=p75;
    end
end

%% sort database in order of descending registration error x 3 and decreasing coverage. 
W = mean([3.*nanmean(meta.rmse)./meta.rmse,...
    meta.gridPointN./nanmean(meta.gridPointN)],2);

W = [double(meta.qc),1./W];

[~,n] = sortrows(W);
meta = structfun(@(x) ( x(n,:) ), meta, 'UniformOutput', false);


%% Add Registered Strips as Anchors
nreg = find(~isnan(meta.rmse));
x = m.x;
y = m.y;
while ~isempty(nreg)
    
    % Check for redundancy    
    if any(N(:))
        
        % subset the DEMs to those that just overlap the rectangle of
        % existing data coverage
        anyN=any(N);
        minNx = x(find(anyN,1,'first'));
        maxNx = x(find(anyN,1,'last'));
        anyN=any(N');
        maxNy=y(find(anyN,1,'first'));
        minNy =y(find(anyN,1,'last'));
        clear anyN
        
        % find files within boundary
        n = find(meta.xmax > minNx & meta.xmin < maxNx & ...
            meta.ymax > minNy & meta.ymin < maxNy);
        
        % loop through files in boundary and find # of new/existing points
        overlap=meta.overlap;
        for i=1:length(n)
            
            % subset of N at data points in this file
            Nsub=N(meta.gridPointInd{n(i)}) > 0;
            
            % count existing data points
            overlap(n(i))=sum(Nsub);
            
            % count new data points
            meta.gridPointN(n(i)) = sum(~Nsub);
            
        end
       
        
        % refresh overlap
        meta.overlap=overlap; overlap=[];
        
        % remove rendundant files (with already 100% coverage)
        redundantFlag = meta.gridPointN < minNewPixels;
        
        if any(redundantFlag)
            
            % append redundant file records to ouput
            [m,meta] = appendRedundant(m,meta,redundantFlag);
            
            % remove redundant files from lists
            meta = structfun(@(x) ( x(~redundantFlag,:) ), meta, 'UniformOutput', false);
            
            fprintf('%d redundant files removed\n',sum(redundantFlag));
           
            nreg = find(~isnan(meta.rmse));
            
            if isempty(nreg); break; end

        end
        
    end

    i = nreg(1);
    
    fprintf('%d anchor files remaining, %d total files remaining \n',length(nreg),length(meta.f));
    fprintf('adding anchor strip: %s, quality=%d, rmse=%.2fm,coverage=%d pixels\n',meta.f{i}, meta.qc(i), meta.rmse(i), meta.gridPointN(i))
        
    mask=cell(1,2);
    if isfield(meta,'maskPolyx')
        if ~isempty(meta.maskPolyx{i})
            mask=[meta.maskPolyx{i}(:) meta.maskPolyy{i}(:)];
        end
    end
            
    N = addStrip2Mosaic( meta.f{i},m,meta.stripDate(i)-dy0,N,...
        meta.dtrans(i,:)',meta.rmse(i),...
        'mergeMethod',mergeMethodReg,...
        'mask',mask,...
        'addErrors',false);
        
    % add this file name
    m.f=[m.f,meta.f(i)];
    m.reg=[m.reg,{'R'}];
    m.N = N;
        
    % add meta data
    m.overlap=[m.overlap,0];
    m.stripDate = [m.stripDate,meta.stripDate(i)];
    m.qcflag = [m.qcflag,meta.qc(i)];
        
        
    if isfield(meta,'maskPolyx')
        m.maskPolyx = [m.maskPolyx,meta.maskPolyx(i)];
        m.maskPolyy = [m.maskPolyy,meta.maskPolyy(i)];
    end
    
    % remove these files from the meta struct
    in = true(size(meta.f)); in(i) = false;
    meta = structfun(@(x) ( x(in,:) ), meta, 'UniformOutput', false);
    
    C(N ~= 0)=int16(1);
    m.C = C;
    m.Nreg = C;
    
    sumN = sum(N(:) > 0);
    
    fprintf('%.2f%% of tile filled\n',sumN./numel(N).*100)
    
    % check if compeletly covered - if so skip
    if numel(N) == sumN
        fprintf('coverage complete, skipping')
        return
    end
    
    nreg = find(~isnan(meta.rmse));
        
end

function [m,meta] = appendRedundant(m,meta,redundantFlag)

% append redundant file records to ouput
m.f=[m.f,meta.f(redundantFlag)'];
m.reg=[m.reg,repmat({'-'},1,sum(redundantFlag))];
m.dtrans=[m.dtrans,nan(3,sum(redundantFlag))];
m.rmse=[m.rmse,nan(1,sum(redundantFlag))];
m.overlap = [m.overlap,meta.overlap(redundantFlag)'];
m.qcflag = [m.qcflag,meta.qc(redundantFlag)'];
if isfield(meta,'maskPolyx')
    m.maskPolyx =  [m.maskPolyx,meta.maskPolyx(redundantFlag)'];
    m.maskPolyy =  [m.maskPolyy,meta.maskPolyy(redundantFlag)'];
end
m.stripDate = [m.stripDate,meta.stripDate(redundantFlag)'];
