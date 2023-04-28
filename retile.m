function B = retile(varargin)
%retile: mosaic a cell array of rectangular tiles into a single array
%B = retile(A)
%B = retile(A, [overlap])
%B = retile(A, [overlap], ['none','linear'])

%input:
%A = cell array
%overlap = a 1 or 2 element array specifying the number of pixels of
%overlap on each inner tile edge. [i=j overlap] or
% [j overlap, i overlap]
% feathering = specify whether or not to feather overlapping sections and
% what method to use (curently, only "linear" is employed). Linear
% feathering will be the default when an overlap is provided.


%Ian Howat, Ohio State University, ihowat@gmail.com

fflag = 0;
Nx = 0;
%% Parse inputs

A = varargin{1};


if nargin == 1;
    A = cell2mat(A);
end

if nargin > 1
    Nx = varargin{2}(1);
    if length(varargin{2}) == 2
        
        Ny = varargin{2}(1);
    else
        Ny = Nx;
    end
    fflag = 1;
end

if nargin > 2
    for i=3:2:nargin
        
        switch lower(varargin{i})
            case 'feathering'
                switch lower(varargin{i+1})
                    case 'none'
                        fflag = 0;
                    case 'linear'
                        fflag = 1;
                end
                %             case 'backgroundvalue'
                %                 backgroundval = varargin{i+1};
        end
    end
end

%% Merge with linear feathering

% merge first row
B = mergecols(A(1,:),Nx,fflag);

sz1 = size(A,1);

if sz1 > 1
    count = 0;
    for i = 2:sz1
        
        if i>1 for k=1:count fprintf('\b'); end; %delete line before
            count = fprintf('Retile: %.2f%s',100*i/sz1,'%');
        end
        
        % merge next row
        C = mergecols(A(i,:),Nx,fflag);
        
        a = B(end-Ny+1:end,:);
        b = C(1:Ny,:);
        
        if fflag
            
            if ~isfloat(a)
                
                a = single(a);
                b = single(b);
                a(a == -9999) = NaN;
                b(b == -9999) = NaN;
            end
            
            % weigting matrix
            I = repmat(linspace(0,1,Ny)',[1,size(B,2)]);
            
            D =  a.*(1-I) + b.*I;
            
            n = isnan(a) & ~isnan(b);
            D(n) = b(n);
            
            n = ~isnan(a) & isnan(b);
            D(n) = a(n);
        else
            
            D =  [a(1:Ny/2,:);b(Ny/2+1:end,:)];
            
        end
        
        if ~isfloat(B)
            
            D(isnan(D)) = -9999;
            D = int16(D);
        end
        % merge row onto array
        B = [B(1:end-Ny,:); D; C(Ny+1:end,:)];
    end
end

disp('done')

function B = mergecols(A,N,fflag)
% merge columns with feathering

B = A{1};

sz2 = size(A,2);

if sz2 > 1
    
    for i = 2:sz2
        
        a = B(:,end-N+1:end);
        b = A{i}(:,1:N);
        
        if fflag
            
            if ~isfloat(a)
                
                a = single(a);
                b = single(b);
                a(a == -9999) = NaN;
                b(b == -9999) = NaN;
            end
            
            % weigting matrix
            I = repmat(linspace(0,1,N),[size(a,1),1]);
            
            D =  a.*(1-I) + b.*I;
            
            n = isnan(a) & ~isnan(b);
            D(n) = b(n);
            
            n = ~isnan(a) & isnan(b);
            D(n) = a(n);
        else
            
            D =  [a(:,1:N/2),b(:,N/2+1:end)];
            
        end
        
        if ~isfloat(B)
            
            D(isnan(D)) = -9999;
            D = int16(D);
            
        end
        % merge row onto array
        B = [B(:,1:end-N), D, A{i}(:,N+1:end)];
    end
    
end
