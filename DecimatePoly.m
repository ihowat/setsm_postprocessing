function [C_out,i_rem]=DecimatePoly(C,opt)
% Reduce the complexity of a 2D simple (i.e. non-self intersecting), closed 
% piecewise linear contour by specifying boundary offset tolerance. 
% IMPORTANT: This function may not preserve the topology of the original 
% polygon. 
%
% INPUT ARGUMENTS:
%   - C     : N-by-2 array of polygon co-ordinates, such that the first, 
%             C(1,:), and last, C(end,:), points are the same. 
%   - opt   : opt can be specified in one of two ways: 
%             ----------------------APPROACH #1 (default) -----------------
%             - opt : opt=[B_tol 1], where B_tol is the maximum acceptible 
%                     offset from the original boundary, B_tol must be 
%                     expressed in the same lenth units as the co-ords in 
%                     C. Default setting is B_tol=Emin/2, where Emin is the
%                     length of the shortest edge.  
%             ----------------------APPROACH #2----------------------------
%              - opt : opt=[P_tol 2], where P_tol is the fraction of the 
%                      total number of polygon's vertices to be retained. 
%                      Accordingly, P_tol must be a real number on the
%                      interval (0,1).
%
% OUTPUT:
%   - C_out : M-by-2 array of polygon coordinates.
%   - i_rem : N-by-1 logical array used to indicate which vertices were
%             removed during decimation.
%
% ALGORITHM:
% 1) For every vertex compute the boundary offset error.  
% 2) Rank all vertics according to the error score from step 1.
% 3) Remove the vertex with the lowest error.
% 4) Recompute and accumulate the errors for the two neighbours adjacent to
%    the deleted vertex and go back to step 2. 
% 5) Repeat step 2 to 4 until no more vertices can be removed or the number
%    of vertices has reached the desired number.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: Jan.2011
%

% Check the input args
if nargin<2, opt={}; end
opt=CheckInputArgs(C,opt);

N=size(C,1);
i_rem=false(N,1); 
if N<=4, 
    C_out=C; 
    return
end

% Tolerance parameter, perimeter and area of the input polygon
[Po,Emin]=PolyPerim(C);
B_tol=Emin/2;
Ao=PolyArea(C);
No=N-1;

if ~isempty(opt), B_tol=opt(1); end

Nmin=3;
if opt(2)==2
    Nmin=round((N-1)*opt(1));
    if (N-1)==Nmin, return; end
    if Nmin<3, Nmin=3; end
end

% Remove the (repeating) end-point
C(end,:)=[];
N=N-1;


% Compute the distance offset errors --------------------------------------
D31=circshift(C,[-1 0])-circshift(C,[1 0]);
D21=C-circshift(C,[1 0]);
dE_new2=sum(D31.^2,2); % length^2 of potential new edges

% Find the closest point to the current vertex on the new edge
t=sum(D21.*D31,2)./dE_new2; 
if t<0, t(t<0)=0; end %#ok<*BDSCI>
if t>1, t(t>1)=1; end
V=circshift(C,[1 0])+bsxfun(@times,t,D31);

% Evaluate the distance^2
Err_D2=sum((V-C).^2,2);

% Initialize distance error accumulation array 
DEAA=zeros(N,1);


% Begin decimation --------------------------------------------------------
idx_ret=1:N; % keep track of retained vertices
while true
    
    % Find the vertices whose removal will satisfy the decimation criterion
    idx_i=Err_D2<B_tol;
    if sum(idx_i)==0 && N>Nmin && opt(2)==2
        B_tol=B_tol*sqrt(1.5);
        continue
    end

    idx_i=find(idx_i);
    if isempty(idx_i) || N==Nmin, break; end
    N=N-1;
    
    % Vertex with the smallest net error
    [~,i_min]=min(Err_D2(idx_i));
    idx_i=idx_i(i_min);

    
    % Update the distance error accumulation array 
    DEAA(idx_i)=DEAA(idx_i)+sqrt(Err_D2(idx_i));
    
    i1=idx_i-1; if i1<1, i1=N; end
    i3=idx_i+1; if i3>N, i3=1; end
    
    DEAA(i1)=DEAA(idx_i);
    DEAA(i3)=DEAA(idx_i);
    
    % Recompute the errors for the vertices neighbouring the vertex marked
    % for deletion
    i1_1=i1-1; if i1_1<1, i1_1=N; end
    i1_3=i3;
    
    i3_1=i1;
    i3_3=i3+1; if i3_3>N, i3_3=1; end
    
    err_D1=RecomputeErrors(C([i1_1,i1,i1_3],:));
    err_D3=RecomputeErrors(C([i3_1,i3,i3_3],:));
    
    % Upadate the errors
    Err_D2(i1)=(sqrt(err_D1)+ DEAA(i1)).^2;
    Err_D2(i3)=(sqrt(err_D3)+ DEAA(i3)).^2;
        
    % Remove the vertex
    C(idx_i,:)=[];
    idx_ret(idx_i)=[];
    DEAA(idx_i)=[];
    
    Err_D2(idx_i)=[];
    
end
C=[C;C(1,:)]; C_out=C;

i_rem(idx_ret)=true;
i_rem=~i_rem;
i_rem(end)=i_rem(1);

% Perimeter and area of the simplified polygon
P=PolyPerim(C);
A=PolyArea(C);

% Performance summary
fprintf('\t\t# of verts\t\tperimeter\t\tarea\n')
fprintf('in\t\t%-5u\t\t\t%-.2f\t\t\t%-.2f\n',No,Po,Ao)
fprintf('out\t\t%-5u\t\t\t%-.2f\t\t\t%-.2f\n',N,P,A)
fprintf('-----------------------------------------------------\n')
fprintf('change\t%-5.2f%%\t\t\t%-5.2f%%\t\t\t%-5.2f%%\n\n',(N-No)/No*100,(P-Po)/Po*100,(A-Ao)/Ao*100)


%==========================================================================
function err_D2=RecomputeErrors(V)
% Recompute the distance offset error for a small subset of polygonal 
% vertices.  
%
%   - V     : 3-by-2 array of triangle vertices, where V(2,:) is the vertex
%             marked for removal.

% Compute the distance offset error ---------------------------------------
D31=V(3,:)-V(1,:);
D21=V(2,:)-V(1,:);
dE_new2=sum(D31.^2,2); % length^2 of potential new edge

% Find the closest point to the current vertex on the new edge
t=sum(D21.*D31,2)/dE_new2; 
if t<0, t(t<0)=0; end
if t>1, t(t>1)=1; end
p=V(1,:)+bsxfun(@times,t,D31);

% Evaluate the distance^2
err_D2=sum((p-V(2,:)).^2);


%==========================================================================
function [P,Emin]=PolyPerim(C)
% Polygon perimeter.

dE=C(2:end,:)-C(1:(end-1),:);
dE=sqrt(sum(dE.^2,2));
P=sum(dE);
Emin=min(dE);


%==========================================================================
function A=PolyArea(C)
% Polygon area. 

dx=C(2:end,1)-C(1:(end-1),1);
dy=C(2:end,2)+C(1:(end-1),2);
A=abs(sum(dx.*dy)/2);


%==========================================================================
function opt=CheckInputArgs(C,opt)
% Check the validity of the input arguments

siz=size(C);
if numel(siz)~=2 || siz(2)>siz(1) || ~isnumeric(C) || ~ismatrix(C) || ndims(C)~=2
    error('First input argument must be a N-by-2 array of polygon vertices')
end

if ~isequal(C(1,:),C(end,:))
    error('First and last points in C must be the same')
end

if isempty(opt), return; end

if ~isnumeric(opt) || numel(opt)~=2
    error('Incorrect entry for 2nd input argument')
end

if ~(opt(2)==1 || opt(2)==2)
    error('Incorrect entry for 2nd input argument. opt(2) must be set to 1 or 2.')
end

if opt(2)==1 && opt(1)<=eps
    error('Incorrect entry for 2nd input argument. When opt(2)==1, opt(1) must be greater than zero.')
end

if opt(2)==2 && (opt(1)<=eps || opt(1)>=1)
    error('Incorrect entry for 2nd input argument. When opt(2)==2, opt(1) must be on the interval (0,1).')
end

