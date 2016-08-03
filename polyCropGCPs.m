function [n]=polyCropGCPs(gcp,xv,yv,varargin)

nonRectPoly=true;
if ~isempty(varargin)
    switch lower(varargin{1})
        case 'rectangle'
            nonRectPoly=false;
    end
end


% fast rectangle crop gcp data to rectangularized DEM boundarys
n= gcp.x >= min(xv) & gcp.x <=max(xv) &...
    gcp.y >= min(yv) & gcp.y <=max(yv);


n=find(n); % reduce index to fast crop selection

if nonRectPoly
    % polygon crop THIS IS NOT NEEDED FOR RECTANGULAR TILES
    nn = inpolygon(gcp.x(n),gcp.y(n),xv,yv); % polygon crop
   
    n=n(nn); clear nn % reduce index to fast crop selection
    
end
