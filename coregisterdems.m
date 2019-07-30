function [z2out,p,perr,d0] = coregisterdems(x1,y1,z1,x2,y2,z2,varargin)
% COREGISTERDEM registers a floating to a reference DEM
%
% [z2r,trans,trans_err,rms] = coregisterdems(x1,y1,z1,x2,y2,z2) registers the
% floating DEM in 2D array z2 with coordinate vectors x2 and y2 to the
% reference DEM in z1 using the iterative procedure in Nuth and Kaab,
% 2011. z2r is the regiestered DEM, p is the [dz,dx,dy] transformation
% parameters, with their 1-sigma errors in trans_err, and rms is the rms of the
% transformation in the vertical from the residuals. If the registration fails
% due to lack of overlap, NaNs are returned in p and perr. If the registration
% fails to converge or exceeds the maximum shift, the median vertical offset is 
% applied.
%
% [...]= coregisterdems(x1,y1,z1,x2,y2,z2,m1,m2) allows a data mask to be applied
% where 0 values will be ignored in the solution.

% Maximum offset allowed
maxp = 15;

x1=x1(:)';
y1=y1(:);

x2=x2(:)';
y2=y2(:);

if length(x1) < 3 || length(y1) < 3 || length(x1) < 3 || length(y1) < 3
    error('minnimum array dimension is 3')
end

interpflag=true;
if (length(x1) == length(x2)) && (length(y1) == length(y2))
	if ~any(x2-x1) && ~any(y2-y1)
	 	interpflag=false;
	end
end

if length(varargin) == 2
    
    m1=varargin{1};
    m2=varargin{2};
    
end

rx = x1(2)-x1(1); % coordinate spacing
p = [0;0;0]; % initial trans variable
pn = p; % iteration variable
perr = p; %regression errors
pnerr = p; %regression errors
d0 = inf; % initial rmse
it = 1; % iteration step

while it

    if interpflag
        % interpolate the floating data to the reference grid
        z2n = interp2(x2 - pn(2),y2 - pn(3),z2 - pn(1),...
            x1,y1,'*linear');
         if exist('m2','var')
            m2n = interp2(x2 - pn(2),y2 - pn(3),single(m2),...
            x1,y1,'*nearest');
            m2n(isnan(m2n)) = 0; % convert back to uint8
            m2n = logical(m2n);
         end
    else
        z2n=z2-pn(1);
        if exist('m2','var'); m2n=m2; end
    end

   interpflag=true;

    % slopes
    [sx,sy] = gradient(z2n,rx);
    sx = -sx;
    
    fprintf('Planimetric Correction Iteration %d ',it)
    
    % difference grids
    dz = z2n - z1;
   
    if exist('m1','var'); dz(~m1) = NaN; end
    if exist('m2','var'); dz(~m2n) = NaN; end
    
    if ~any(~isnan(dz(:))); disp('No overlap'); z2out=z2; p=[NaN;NaN;NaN]; d0=NaN; return; end
	 
    % filter NaNs and outliers
    n = ~isnan(sx) & ~isnan(sy) & ...
        abs(dz - nanmedian(dz(:))) <= 3*nanstd(dz(:));
    
    if sum(n(:)) < 10 
        fprintf('Too few (%d) registration points, quitting and returning NaNs\n',sum(n(:))); 
       	z2out=z2; p=[NaN;NaN;NaN]; perr=[NaN;NaN;NaN]; d0=[NaN];  
        return 
    end
    
    % get RMSE and break if below threshold
    d1 = sqrt(mean(dz(n).^2));
    
    % keep median dz if first iteration
    if it == 1
        meddz=median(dz(n));
	meddz_err=std(dz(n))./sqrt(sum(n)); 
        d00=sqrt(mean((dz(n)-meddz).^2));
    end
    
    fprintf('rmse= %.3f ',d1)
    
    if d0 - d1 < .001 || isnan(d0) || it == 5
        
        fprintf('stopping \n')
        % if fails after first registration attempt, set dx and dy to zero
        % and subtract the median offset
        if it == 2
            fprintf('regression failure, returning median vertical offset: %.3f\n',meddz)
            p(1)=meddz; d0=d00;
            z2out = z2-meddz; 
	elseif it == 5
            fprintf('maximum number of iterations reached\n')
        else
            fprintf('improvement minimum reached\n')
        end
        break
    end
        
    %keep this adjustment
    p = pn;
    perr = pnerr;
    d0 = d1;
    z2out = z2n;
   
    % build design matrix
    X = [ones(size(dz(n))),sx(n),sy(n)];

    % solve for new adustment
    p1 = X\dz(n);
    
    % calculate p errors
    [~,R,perm] = qr(X,0);
    RI = R\eye(3);
    nu = size(X,1)-size(X,2); % Residual degrees of freedom
    yhat = X*p1;              % Predicted responses at each data point.
    r = dz(n)-yhat;           % Residuals.
    normr = norm(r);
    
    rmse = normr/sqrt(nu);      % Root mean square error.
    tval = tinv((1-0.32/2),nu);
    
    se = zeros(size(X,2),1);
    se(perm,:) = rmse*sqrt(sum(abs(RI).^2,2));
    p1err = tval*se;
    
    % update shifts
    pn = p + p1;
    pnerr = sqrt(perr.^2 + p1err.^2);
    
    % display offsets
    fprintf('offset(z,x,y): %.3f, %.3f, %.3f\n',pn')
    
    if any(abs(pn(2:end)) > maxp)
        fprintf('maximum horizontal offset reached, returning median vertical offset: %.3f\n',meddz)
        p=[meddz;0;0]; perr=[meddz_err;0;0]; d0=d00;
        z2out = z2-meddz;
        break
    end

    % update iteration vars
    it = it+1;
end


