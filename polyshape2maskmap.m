function [BW,x,y]=polyshape2maskmap(varargin)
% POLYSHAPE2MASKMAP convert a polyshape to a bit mask in map coordinates
%
% BW = polyshape2maskmap(p,x,y) where p is the polyshape with vertices in
% the same coordinates system map coordinate vectors x and y. BW is the 
% length(y) by length(x) logical array with true for regions and false for
% backgroun or holes.
%
% BW = polyshape2maskmap(p,[x0 x1 y0 y1],[dx dy]) uses the ranges x0 to x1
% and y0 to y1 and the pixel spacing dx and dy to build the coordinate
% vectors and array.
% 
% [BW,x,y] = polyshape2maskmap(...) returns the x and y coordinate vectors

% Ian Howat '06-Aug-2020 11:38:23'

p = varargin{1};

% if length of second input is 4, assume its ranges and make map
% coordinates
if length(varargin{2}) == 4

    if length(varargin{3}) == 1
        varargin{3}(2) = varargin{3}(1);
    end
    
    x = varargin{2}(1): varargin{3}(1):varargin{2}(2);
    y = varargin{2}(4):-varargin{3}(2):varargin{2}(3);
    
else
    x = varargin{2};
    y = varargin{3};
end
    
% make bitmask grid
BW = false(length(y),length(x));

if isempty(p.Vertices)
%    error('polyshape has no vertices')
    fprintf('polyshape has no vertices\n')
    return
end

% index the polygons within the vertices vectors
NR = [0;find(isnan(p.Vertices(:,1)));length(p.Vertices(:,1))+1];

% locate the holes
TF = ishole(p);

% loop through polygons
for nr = 1:length(NR)-1
    
    % pull out this polygon's vertices
    xp = p.Vertices(NR(nr)+1:NR(nr+1)-1,1);
    yp = p.Vertices(NR(nr)+1:NR(nr+1)-1,2);
   
    if TF(nr) %its a hole
        % set BW within polygon to false
        BW(roipoly(x,y,BW,xp,yp))=false;
    else %its a region
         % set BW within polygon to true
        BW(roipoly(x,y,BW,xp,yp))=true;
    end
end
