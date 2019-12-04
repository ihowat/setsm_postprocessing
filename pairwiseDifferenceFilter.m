function [f,P] = pairwiseDifferenceFilter(z,varargin)

sz = size(z);

prctileThresh = 25;
epsilon=10;
minpts=3;

N = sum(~isnan(z),3);
mask = N > minpts;

if length(varargin) > 1
    
    varargin(cellfun(@isstr,varargin))=...
        lower(varargin(cellfun(@isstr,varargin)));
    
    narg=find(strcmp('mask',varargin));
    if narg
        mask=varargin{narg+1};
        
        mask = mask & N > minpts;
        fprintf('Applying mask\n')
        
    end
    
    narg=find(strcmp('minn',varargin));
    if narg
        minN=varargin{narg+1};
        
        % return number of repeats/pixel
        fprintf('Applying minN = %.1f\n',minN)
        
        mask = mask & N > minN;
        
    end
    
    narg=find(strcmp('minmad',varargin));
    if narg
        minMad=varargin{narg+1};
        
        fprintf('Applying minMad = %.1f\n',minMad)
        % calculate MAD of each pixel
        z_mad = mad(z,1,3);
        
        mask = mask & z_mad > minMad;
    end
    
    narg=find(strcmp('epsilon',varargin));
    if narg
        epsilon=varargin{narg+1};
    end
    
    narg=find(strcmp('minpts',varargin));
    if narg
        minpts=varargin{narg+1};
    end
    
    narg=find(strcmp('prctilethresh',varargin));
    if narg
        prctileThresh=varargin{narg+1};
    end
    
    narg=find(strcmp('datenum',varargin));
    if narg
        t=varargin{narg+1};
        p = nan(size(z,1),size(z,2),2);
        
    end
    
end

fprintf('Epsilon = %.1f\n',epsilon)
fprintf('minpts = %d\n',minpts)


[j,k]  = find(mask);

sz = size(z);
f = true(sz);
Npix = length(j);
i=1;

fprintf('%d pixles\n',Npix)

if ~isinf(prctileThresh)
    
    fprintf('clustering pairwise differences, prctileThresh = %d\n',prctileThresh)
    
    for i=1:Npix
        
        zp = reshape(z(j(i),k(i),:),sz(3),1,1);
        
        n = find(~isnan(zp));
        
        zp = zp(n);
        
        if exist('t','var')
            p0 = polyfit(t(n),zp,1);
            zp = zp - polyval(p0,t(n));
            P(j(i),k(i),1) = p0(1);
             P(j(i),k(i),2) = p0(2);
        end
        
        dzp = repmat(zp',length(n),1) - repmat(zp,1,length(n));
        
        dzp(diag(true(length(n),1))) = NaN;
        
        dzp_abs_25pt = prctile(abs(dzp),prctileThresh)';
        
        idx = dbscan([dzp_abs_25pt,zp],epsilon,minpts);
        
        [~,minN] = min(dzp_abs_25pt);
        
        idx = idx == idx(minN);
        
        n = n(~idx);
        
        f(j(i),k(i),n) = false;
        
    end
    
else
    
    fprintf('clustering point values\n')
    for i=1:Npix
        zp = reshape(z(j(i),k(i),:),sz(3),1,1);
        
        n = find(~isnan(zp));
        
        zp = zp(n);
        
        if exist('t','var')
            p = polyfit(t(n),zp,1);
            zp = zp - polyval(p,t(n));
        end

        idx = dbscan(zp,epsilon,minpts);
        
        if ~any(idx ~= -1)
            
            %             idx = dbscan(zp,2.*epsilon,minpts);
            %
            %             if ~any(idx ~= -1)
            
            continue
        end
        
    end
    
    n = n(idx == -1);
    
    
    f(j(i),k(i),n) = false;
    
end