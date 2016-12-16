function batchAlignTile(f,freg)
% batchAlignTile: recursively apply alignTile to all nregistered tiles
%
% batchAlignTile(f,freg) where f and freg are cellstr of unregistered and
% registered files, respectively. Filenames must be in the format
% /rr_cc_*.mat where rr and cc are the tile row and column numbers.

% orient horizontally
f=f(:)';
freg=freg(:)';

while ~isempty(f)
    
    fprintf('%d unregistered tiles, %d registered tiles\n',length(f),length(freg));
    
    % get file names for row/col extract
    [~,fname]=cellfun(@fileparts, f, 'UniformOutput',false);
    [~,fregname]=cellfun(@fileparts, freg, 'UniformOutput',false);
    
    % make sure no reg'd unregs
    [~,ia]=intersect(fname,fregname);
    fname(ia)=[];
    f(ia)=[];
    
    % tile row/col index from filename
    tc = char(fname(:));
    tr= str2double(tc(:,1:2));
    tc= str2double(tc(:,4:5));
    
    tcreg = char(fregname(:));
    trreg= str2double(tcreg(:,1:2));
    tcreg= str2double(tcreg(:,4:5));
    
    % get up/down/right/left neighbors
    A = [[0 1];[1 0];[-1 0];[0 -1]];
    
    % for each unregistered file, get neighboring registered files
    nreg=zeros(length(f),4);
    for i = 1:length(f)
        
        for j=1:4
            
            n = find(tr(i)+A(j,1) == trreg & tc(i)+A(j,2) == tcreg);
            if ~isempty(n)
                nreg(i,j) = n;
            end
            
        end
    end
    
    % count number of registered files
    Nreg=sum(nreg ~= 0,2);
    
    % make surge there are neighbors
    if ~any(Nreg)
        fprintf('%d unregistered tiles remain with no registered neigbors\n',length(f));
        break
    end
    
    % select tile with most registered  neighbors (or first in list of tie).
    [~,n]=max(Nreg);
    n=n(1);
    
    nreg=nreg(n,nreg(n,:) ~= 0);
    fprintf('registering %s to %s, %s, %s, %s',f{n},freg{nreg}); fprintf('\n');
    
    [outFlag,outname] = alignTile({f{n},freg{nreg}});
    
    if outFlag
        f(n)=[];
        freg=[freg {outname}];
        failCount = 0;
    else
        failCount=failCount + 1;
        if failCount == length(f)
            fprintf('%d tiles cannot be registered, quitting\n',length(f));
            break
        end
            
    end
    
end








    
