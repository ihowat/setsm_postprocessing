function A = tile(B,nx,ny,varargin)
%tile: split a 2D array into equally-sized subarrays.
%A = tile(B,nx,ny,overlap)
%input:
%B = 2D array
%nx = number of tiles in column direction
%ny = number of columns in row direction
%overlap = ammount of overlap of tiles in pixels
%OUTPUT
%A = cell array of size ny by nx with each cell containing a tile.

%Ian Howat, Ohio State University, ihowat@gmail.com

%set overlaps defaults
overlapx = 0;
overlapy = 0;

%parse inputs
switch length(varargin)
    case 1
        %1 overlap given - set to both
        overlapx = varargin{1};
        overlapy = varargin{1};
    case 2
        % x and y overlaps given
        overlapx = varargin{1};
        overlapy = varargin{2};
end

% array dimensions
nl = size(B,1);
ns = size(B,2);


if ny > 1
 l = floor(linspace(1,nl,ny+1));
 minl = [1,l(2:end-1)+1];
 maxl = [l(2:end-1),nl];
 dl = diff([minl;maxl]);
 minl(2:end) = minl(2:end) - overlapy/2; 
 maxl(1:end-1) = maxl(1:end-1) + overlapy/2;
 minl(minl < 1) = 1;
 maxl(maxl > nl) = nl;
else
   maxl = nl;
   minl = 1;
end

if nx > 1
s = floor(linspace(1,ns,nx+1));
mins = [1,s(2:end-1)+1];
maxs = [s(2:end-1),ns];
ds = diff([mins;maxs]);
mins(2:end) = mins(2:end) - overlapx/2;
maxs(1:end-1) = maxs(1:end-1) + overlapx/2;
mins(mins < 1) = 1;
maxs(maxs > ns) = ns;
else 
    maxs = ns;
    mins = 1;
end

A = cell(ny,nx);

% subset loop
for i=1:ny
   for j=1:nx
   
    A{i,j} =  B(minl(i):maxl(i),mins(j):maxs(j));

   end
end


