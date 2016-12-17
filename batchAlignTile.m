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
    
    fprintf('%d unregistered tiles, %d registered tiles, %d failed to register\n',length(f),length(freg),failCount);
    
    % get file names for row/col extract
    [~,fname]=cellfun(@fileparts, f, 'UniformOutput',false);
    [~,fregname]=cellfun(@fileparts, freg, 'UniformOutput',false);
    
    % make sure no reg'd unregs
    [~,ia]=intersect(fname,fregname);
    fname(ia)=[];
    f(ia)=[];
    
    % tile row/col index from filename
    tc = char(fname(:));
    tr= str2num(tc(:,1:2));
    tc= str2num(tc(:,4:5));
    
    tcreg = char(fregname(:));
    trreg= str2num(tcreg(:,1:2));
    tcreg= str2num(tcreg(:,4:5));
    
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
    
    % rank tiles in descending order of number of neighboring registered tiles
    [~,n]=sort(Nreg,'descend');
    
    % initialize fail counter
    failCount=0
    
    % sequentially attempt each file by n, breaking if successful or all fail.
    while failCount < length(f)
    
        % extract this n
        n1=n(1 + failCount);
        nreg1=nreg(n,nreg(n,:) ~= 0);
    
        fprintf('registering %s to %s, %s, %s, %s',f{n1},freg{nreg1}); fprintf('\n');
        
        % try to coregister
        [outFlag,outname] = alignTile({f{n1},freg{nreg1}});
    
        % outFlag = true is success, so update lists, break and redo neignbor search
        if outFlag
            f(n1)=[]; % remove this file from the unreg list 
            freg=[freg {outname}]; % add this file to the reg list
            break
        else % if coreg failed update count and try the next one in the n list
            failCount=failCount + 1; 
            fprintf('%d failed tiles, %d unregistered tiles left\n',failCount, length(f))
        end  
    end
    
    % if all remaining unreg tiles failed, quit
    if failCount >= length(f)
        fprintf('%d tiles cannot be registered, quitting\n',length(f));
        break
    end
end








    
