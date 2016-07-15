function f1 = orderPairs(demdir,f0)

% original pair index
n0 = 1:length(f0);

% pair boundary vertices
xv = zeros(5,length(n0));
yv = zeros(5,length(n0));

% imrect format vertices
R0 = zeros(length(n0),4);

% ordered pair index
n1 = zeros(1,length(n0));

% loop through files and get border information
for i=1:length(f0);
    try
	d = readGeotiff([demdir,'/',f0{i}],'mapinfoonly');
    catch 
	error(['error reading map info from ',demdir,'/',f0{i}])
    end
    R0(i,:) = [min(d.x),min(d.y),max(d.x)-min(d.x),max(d.y)-min(d.y)];
    xv(:,i) = [min(d.x); min(d.x); max(d.x); max(d.x); min(d.x)];
    yv(:,i) = [min(d.y); max(d.y); max(d.y); min(d.y); min(d.y)];
    clear d
end

% get aspect ratio of strip
ar = (max(R0(:,1) + R0(:,3)) - min(R0(:,1)))./...
    (max(R0(:,2) + R0(:,4)) - min(R0(:,2)));

% start with bottom left pair
if ar >= 1;
    [~,n1(1)] = min(R0(:,1));
else
    [~,n1(1)] = min(R0(:,2));
end

% put this starting pair into the ordered index
R1  = R0(n1(1),:);

n0(n1(1)) =[];
R0(n1(1),:) =[];

% loop through pairs and sequentially add the next pair with the most
% overlap
i=2;
for i=2:length(f0);
    
    % loop through all scenes and calculate overlap
    A=zeros(length(n0),1);
    j=1; for j=1:length(A); A(j) = rectint(R1,R0(j,:)); end
    
    if ~any(A > 0);
        disp('Break in overlap detected, returning this segment only')
        break
    end
    
    % find pair with max overlap to pair i
    [~,mxAn] = max(A);
    n1(i) = n0(mxAn);
    
    % expand the combined footprint to add new pair
    R0xmax = R0(mxAn,1)+ R0(mxAn,3);
    R0ymax = R0(mxAn,2)+ R0(mxAn,4);
    
    R1xmax = R1(1)+ R1(3);
    R1ymax = R1(2)+ R1(4);
    
    R1(1) = min([R1(1),R0(mxAn,1)]);
    R1(2) = min([R1(2),R0(mxAn,2)]);
    R1(3) = max([R1xmax,R0xmax])-R1(1);
    R1(4) = max([R1ymax,R0ymax])-R1(2);
    
    % remove this pair from the list
    n0(mxAn) =[];
    R0(mxAn,:) =[];
end

% remove files with no overlap
n1(n1==0) = [];
f1 = f0(n1);



