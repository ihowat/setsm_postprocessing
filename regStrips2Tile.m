function [meta,m,N,C] = regStrips2Tile(meta,m,N,C,mergeMethodReg,dy0)
% regStrips2Tile add registered strips to tile
%
% [meta,m,N,C] = regStrips2Tile(meta,m,N,C,,mergeMethodReg) 
% determines if registration data exist for strips in meta.f, either in the 
% metastruct (fields dtrans and rmse) or in *_reg.txt files, and passes the
% the strips to the mosaicker. dy0 is the zero day for the image date
% (ArcticDEM uses Y2K)
%
% Ian Howat, ihowat@gmail.com, Ohio State

% if dtrans is not in the meta, look for reg files.
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

%filter out data over 4m std
meta.dtrans(meta.rmse > 4,:) = NaN;
meta.rmse(meta.rmse > 4) = NaN;

% sort database by p75
[~,n] = sort(meta.rmse,'ascend');
meta = structfun(@(x) ( x(n,:) ), meta, 'UniformOutput', false);

nreg=find(~isnan(meta.rmse));

%% Add Registered Strips as Anchors
if ~isempty(nreg)
    
    for i=1:length(nreg)

        fprintf('adding anchor strip %d of %d: %s\n',i,length(nreg),meta.f{i})
        if isfield(meta,'maskPolyx')
            mask=[meta.maskPolyx{i}(:) meta.maskPolyy{i}(:)];
        else
            mask=cell(1,2);
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
        m.maskPolyx = [m.maskPolyx,meta.maskPolyx(i)];
        m.maskPolyy = [m.maskPolyy,meta.maskPolyy(i)];
        
    end
    
    % remove these files from the meta struct
    in=length(nreg)+1:length(meta.f);
    meta = structfun(@(x) ( x(in,:) ), meta, 'UniformOutput', false);
    
    C(N ~= 0)=int16(1);
    m.C = C;
    m.Nreg = C;
    
end
