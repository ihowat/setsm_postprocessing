function dZ = getOffsets(subTileFiles,subTileRow,subTileCol,buff)
% getOffsets returns a set of vertical adjustments based on tile buffers
%
% dZ = getOffsets(subTileFiles,subTileRow,subTileCol,buff) returns a
% vector dZ of vertical offsets for the cell list of subtile file names
% subTileFiles determined by weighted LSQ adjustment of the differences
% over subtile buffers  of  width buff.  The subTileRow and subTileCol
% vectors determine neighbors.

NsubTileFiles = length(subTileFiles);

% output file for coregistration offsets so that we don't need to
% recalculate them all if we want to change something
% assume will be insignificant difference between 10m and 2m adjusments, so
% can reuse
subTileDir=fileparts(subTileFiles{1});
regFile=[subTileDir,'/tileReg.mat'];

%This would make resolution-dependent tile adjustment
%[subTileDir,subTileName]=fileparts(subTileFiles{1});
% resStr=strsplit(subTileName,'_');
% resStr=resStr{end};
%regFile=[subTileDir,'/tileReg',resStr,'.mat'];

% check if coregistration file exists and load it if so
if exist(regFile,'file')
    load(regFile)
else
    % no coregistration file, so calculate all subtile offsets
    
    % iteration output variables
    nrt = nan(NsubTileFiles-1,1);
    dzup = nan(NsubTileFiles-1,1);
    dzrt = nan(NsubTileFiles-1,1);
    dzup_mad = nan(NsubTileFiles-1,1);
    dzrt_mad = nan(NsubTileFiles-1,1);
    N=nan(NsubTileFiles,1);
    
    for n = 1:NsubTileFiles-1
        
        fprintf('subtile %d of %d ',n,NsubTileFiles-1)
        
        % make sure this subtile has a za_med var
        if  ~ismember('za_med',who('-file',subTileFiles{n}))
            fprintf('variable za_med doesn''t exist, skipping\n')
            continue
        end
        
        m0 = matfile(subTileFiles{n});
        N(n) = max(max(m0.N));
        
        
        thisRow=subTileRow(n);
        thisCol=subTileCol(n);
     
        % index of subtile above
        n1 = find(subTileRow == thisRow+1 & subTileCol == thisCol);
        
        if ~isempty(n1)
                 
            % check up neighbor has za_med var
            if ismember('za_med',who('-file',subTileFiles{n1}))
                
                % load top buffer of the bottom (nth) subtile of pair
                % check to make sure buffer is > 10% land
                l0 = m0.land(2:2*buff,2:end-1);
                if sum(l0(:))/numel(l0) > 0.1
                    z0 = m0.za_med(2:2*buff,2:end-1);
                    z0(~l0) = NaN;
                    y0 = m0.y(2:2*buff,:);
                    x0 = m0.x(:,2:end-1);
                    
                    % load bottom buffer of the top (nth+1) subtile of pair
                    m1 = matfile(subTileFiles{n1});
                    z1 = m1.za_med(end-(2*buff-1):end-1,2:end-1);
                    y1 = m1.y(end-(2*buff-1):end-1,:);
                    x1 = m1.x(:,2:end-1);
                    
                    % check dimensions of z0 and z0 buffers consistent
                    if ~any(y0 ~= y1) && ~any(x0 ~= x1)
                        dzn = z0(:)-z1(:);
                        if sum(~isnan(dzn))/numel(dzn) > 0.1
                            dzup(n) = nanmedian(dzn);
                            dzup_mad(n) = mad(dzn,1);
                        end
                    else
                        warning('inconsistent buffers, offset not computed')
                    end
                end
            end
        end
        
        % check right neighbor exists
       % nrtn = find(subTileNum(n)+Nrows ==  subTileNum);
       nrtn = find(subTileRow == thisRow & subTileCol == thisCol+1);
        
        if isempty(nrtn)
            fprintf('\n')
            continue
        end
        
        nrt(n) = nrtn;
        
        % check right neighbor has za_med var
        if ismember('za_med',who('-file',subTileFiles{nrt(n)}))
            
            % check to make sure buffer is > 10% land
            l0 = m0.land(2:end-1,end-(2*buff-1):end-1);
            if sum(l0(:))/numel(l0) > 0.1
                
                % load right buffer of the left (nth) subtile of pair
                z0 = m0.za_med(2:end-1,end-(2*buff-1):end-1);
                z0(~l0) = NaN;
                y0 = m0.y(2:end-1,:);
                x0 = m0.x(:,end-(2*buff-1):end-1);
                
                % load left buffer of the right subtile of pair
                m1 = matfile(subTileFiles{nrt(n)});
                z1 = m1.za_med(2:end-1,2:(2*buff));
                y1 = m1.y(2:end-1,:);
                x1 = m1.x(:,2:(2*buff));
                
                % check dimensions of z0 and z0 buffers consistent
                if ~any(y0 ~= y1) && ~any(x0 ~= x1)
                    dzn = z0(:)-z1(:);
                    if sum(~isnan(dzn))/numel(dzn) > 0.1
                        dzrt(n) = nanmedian(dzn);
                        dzrt_mad(n) = mad(dzn,1);
                    end
                else
                    warning('inconsistent buffers, offset not computed')
                end
            end
        end
        
        fprintf('\n')
    end
    
    % get N for the last file
    n = NsubTileFiles;
    if ismember('N',who('-file',subTileFiles{n}))
        m0 = matfile(subTileFiles{n});
        N(n) = max(max(m0.N));
    end
    
    save(regFile,'nrt','dzup','dzrt','dzup_mad','dzrt_mad','N')
    
end

% adjustment
%build pair indexes
n1 = [(1:NsubTileFiles-1)';(1:NsubTileFiles-1)'];
n2 = [(2:NsubTileFiles)';nrt];

% concoctenate offsets/errors
dz = [dzup;dzrt];
dze = [dzup_mad;dzrt_mad];

% remove nan pairs or with large offsets and/or mad values
n = any(isnan([n1 n2 dz dze]),2);
n = n | abs(dz) > 50 | dze > 2;

n1(n) =[];
n2(n) =[];
dz(n) =[];
dze(n) =[];

Npairs = length(n1);

% Build design and weight matrices
A = zeros(Npairs,NsubTileFiles); % initialize design matrix

linearInd = sub2ind([Npairs NsubTileFiles], (1:Npairs)', n1);
A(linearInd) = 1;
linearInd = sub2ind([Npairs NsubTileFiles], (1:Npairs)', n2);
A(linearInd) = -1;

% locate missing tiles
n_missing  = ~any(A) | isnan(N)';
%
% remove missing tiles
A(:,n_missing) = [];

% % add delta=0 lines
A = [A;diag(ones(1,size(A,2)))];
dz = [dz;zeros(size(A,2),1)];
% dze = [dze;ones(size(A,2),1).*4];

dze = [ones(size(dze));ones(size(A,2),1).*100];

% dze = ones(size(dze),'single');
% dze = [dze;1./sqrt(N(~n_missing))];

% calculate weights
wz = 1./dze.^2;

% build offdset vector
dZ = nan(NsubTileFiles,1);

fprintf('performing LSQ adjustment\n')

dZ(~n_missing) = (wz.*A)\(wz.*dz);

save(regFile,'dZ','-append');
