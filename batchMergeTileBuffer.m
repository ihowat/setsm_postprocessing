function batchMergeTileBuffer(f)
% batchMergeTileBuffer: mergeTileBuffer to all neighboring files
%
% batchMergeTileBuffer(f) where f is a cellstr of files to be merged.
%   Filenames must be in the format /rr_cc_*.mat where rr and cc are the 
%   tile row and column numbers.


[~,fname]=cellfun(@fileparts, f, 'UniformOutput',false);

% tile row/col index from filename
tc = char(fname(:));
tr= str2num(tc(:,1:2));
tc= str2num(tc(:,4:5));

% get up/down/right/left neighbors
A = [[0 1];[1 0];[-1 0];[0 -1]];

% pair counter
c=0;

% pair index variables
n0=0;
n1=0;

% loop through each file and see if neighbors exist, only keeping unique
% pairs
i=1;
for i=1:length(f);
    
    % this row_col string
    si=[num2str(tr(i),'%02d'),'_',num2str(tc(i),'%02d')];
    
    j=1;
    for j=1:4
        
        % that row_col string
        sj=[num2str(tr(i)+A(j,1),'%02d'),'_',num2str(tc(i)+A(j,2),'%02d')];
        
        % make filename
        fj=strrep(f{i},si,sj);
        
        % test if it exists in list
        n = find(strcmp(fj,f));
        
        if ~isempty(n)
            
            % exists, but is it aleady matched?
            [~,ia]=intersect([n i],[n0 n1],'rows');
            
            if ~isempty(ia); continue; end
            
            % new pair, add to count and vectors
            
            c=c+1;
            n0(c,1)=i;
            n1(c,1)=n;
            
        end
        
    end
end

% run mergeTileBuffer on each pair in list
for i=1:length(n0); mergeTileBuffer(f{n0(i)},f{n1(i)}); end
    
    
    




