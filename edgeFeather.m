function W = edgeFeather(M0,M1,varargin)
% EDGEFEATHER get weights for merging two images with edge feathering
%
%W = edgeFeather(M0,M1,buff) where M0 and M1 are equally sized BW images
%with ones for data and zeros for no data. buff is the cutline buffer for
%merging - greater == means more overlap used for feather

buff=0;
f0=0; % overlap fraction where back z weight goes to zero
f1=1; % overlap fraction where back z weight goes to one

maxInterpPixels = 1e6;

% parse input args
for i=1:2:length(varargin)
    
    switch lower(varargin{i})
        
        case 'buffer'
            
            buff=varargin{i+1};
            
        case 'overlaprange'
            
            f0=varargin{i+1}(1);
            f1=varargin{i+1}(2);
                 
        otherwise
            
            error('Unknown input arguments')
    end
end

% crop to just region of overlap
A= single(M0 & M1);

% if no overlap, send 1/0 mask back
if ~any(A(:))
        W=nan(size(M0));
        W(M0) = 1;
        W(M1) = 0;
        return
end

% crop to area of overlap
A(A==0)=NaN;
[~,rb,cb] = cropnans(A,buff);

%Make overlap mask removing isolated pixels
A=single(bwareaopen(...
    ~M0(rb(1):rb(2),cb(1):cb(2)) & M1(rb(1):rb(2),cb(1):cb(2)),...
    1000)); % nodata in back and data in front = 1

A(bwareaopen(...
    M0(rb(1):rb(2),cb(1):cb(2)) & ~M1(rb(1):rb(2),cb(1):cb(2)),...
    1000)) =2 ;  % data in back and nodata in front = 2


% do weight interpolation on low res grid for speed
if numel(A) <= maxInterpPixels; maxInterpPixels=numel(A); end
Ar=imresize(A,maxInterpPixels./numel(A),'nearest');
%Ar=imresize(A,.1,'nearest');

fprintf('%d pixels\n',numel(Ar));

% interpolation grid
[C,R] = meshgrid(1:size(Ar,2),1:size(Ar,1));

% pixles on outside of boundary of overlap region
B = bwboundaries(Ar~=0, 8, 'noholes');
B = cell2mat(B);

if size(B,2) ~= 2
    warning('no overlap boundaries returned, no feathering applied');
    W=nan(size(M0));
    W(M1) = 0;
    W(M0) = 1;
    return
end

n=sub2ind(size(Ar),B(:,1),B(:,2));

warning off
F = scatteredInterpolant(C(n),R(n),double(Ar(n)));
warning on

try
    Ar(Ar==0)=F(C(Ar==0),R(Ar==0));
catch
    warning('interpolation failed, no feathering applied');
    W=nan(size(M0));
    W(M1) = 0;
    W(M0) = 1;
    return
end

Ar=imresize(Ar,size(A),'bilinear');
Ar(A==1 & Ar ~=1)=1;
Ar(A==2 & Ar ~=2)=2;
A=Ar-1;
A(A < 0) = 0;
A(A > 1) = 1;

W=single(M0);
W(rb(1):rb(2),cb(1):cb(2)) = A;
W(M0 & ~M1) = NaN;

clear A

% shift weights so that more of the reference layer is kept
W=(1/(f1-f0)).*W-f0/(f1-f0);
W(W > 1) = 1;
W(W < 0) = 0;

function [A,r,c] = cropnans(varargin)
% cropnans crop array of bordering nans
%
% [A,r,c] = cropnans(A)

A=varargin{1};
buff=0;
if nargin == 2; buff=varargin{2}; end

r = [];
c = [];

M = ~isnan(A);

if ~any(M(:)); return; end

rowsum = sum(M) ~= 0;
colsum = sum(M,2) ~= 0;

c(1) = find(rowsum,1,'first')-buff;
c(2) = find(rowsum,1,'last')+buff;

r(1) = find(colsum,1,'first')-buff;
r(2) = find(colsum,1,'last')+buff;

if c(1) < 1; c(1)=1; end
if r(1) < 1; r(1)=1; end
if c(2) > size(A,2); c(2)=size(A,2); end
if r(2) > size(A,1); r(2)=size(A,1); end

A = A(r(1):r(2),c(1):c(2));
