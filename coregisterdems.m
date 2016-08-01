function [z2out,p,d0] = coregisterdems(x1,y1,z1,x2,y2,z2,varargin)
% COREGISTERDEM registers a floating to a reference DEM
%
% [z2r,trans,rms] = coregisterdems(x1,y1,z1,x2,y2,z2) registers the
% floating DEM in 2D array z2 with coordinate vectors x2 and y2 to the
% reference DEM in z1 using the iterative procedure in Nuth and Kaab,
% 2010. z2r is the regiestered DEM, p is the z,x,y transformation
% parameters and rms is the rms of the transformation in the vertical. 

x1=x1(:)';
y1=y1(:);

x2=x2(:)';
y2=y2(:);


interpflag=true;
if (length(x1) == length(x2)) && (length(y1) == length(y2))
	if ~any(x2-x1) && ~any(y2-y1)
	 	interpflag=false;
	end
end

if length(varargin) == 2;
    
    m1=varargin{1};
    m2=varargin{2};
    
end


rx = x1(2)-x1(1); % coordinate spacing
p = [0;0;0]; % initial trans variable
pn = p; % iteration variable
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
   
    if exist('m1','var') && exist('m2','var'); dz(~m2n | ~m1) = NaN; end
    
    if ~any(~isnan(dz(:))); disp('No overlap'); z2out=z2; p=[NaN;NaN;NaN]; d0=[NaN]; return; end
	 
    % filter NaNs and outliers
    n = ~isnan(sx) & ~isnan(sy) & ...
        abs(dz - nanmedian(dz(:))) <= nanstd(dz(:));
    
    % get RMSE and break if below threshold
    d1 = sqrt(mean(dz(n).^2));
    
    % keep median dz if first iteration
    if it == 1; 
        meddz=median(dz(n)); 
        d00=sqrt(mean((dz(n)-meddz).^2));
    end
    
    fprintf('rmse= %.3f ',d1)
    
    if d0 - d1 < .001 || isnan(d0)
        
        fprintf('stopping \n')
        % if fails after first registration attempt, set dx and dy to zero
        % and subtract the median offset
        if it == 2
            fprintf('regression failure, returning median vertical offset: %.3f\n',meddz)
            p(1)=meddz; d0=d00;
            z2out = z2-meddz; 
        end
        break
    end
    
    %keep this adjustment
    p = pn;
    d0 = d1;
    z2out = z2n;
   
    % build design matrix
    X = [ones(size(dz(n))),sx(n),sy(n)];

    % solve for new adustment
    pn = p + X\dz(n);

    % display offsets
    fprintf('offset(z,x,y): %.3f, %.3f, %.3f\n',pn')
    
    % update iteration vars
    it = it+1;
end

    
% mask = interp2(x2 - pn(2),y2 - pn(3),double(mask),...
%         x1(:)',y1(:),'*linear');
% mask = mask >= 0.1;
    

%clear sx sy x2 y2  dz n d1 d0 it X p z2n maskn
% 
% %% x,y,z - dependent bias
% disp('Ramp and Elevation Bias Correction')
% 
% % grid coordinates for solver
% [X1,Y1] = meshgrid(x1,y1);
% 
% % difference grids
% dz = z2 - z1;
% 
% % filter NaNs for speed and apply std threshold
%  n = abs(dz - nanmedian(dz(:))) <= nanstd(dz(:)) & mask;
% 
% % build design matrix
% X = [ones(size(dz(n))),z1(n),X1(n),Y1(n)];
% 
% % fit
% p = X\dz(n);
% 
% % display offsets
% disp(['y0,dz(x,y,z):',num2str(p')])
% 
% % apply model
% z2 = z2 - (p(1) + p(2).*z2 + p(3).*X1 + p(4).*Y1);
% 
% % recalc differences
% dz = z2 - z1;
% 
% % display final rmse
% rmsfinal = rms(dz(~isnan(dz)));
% fprintf('rmse = %.2f\n',rmsfinal)
% 
% pelev = p;




