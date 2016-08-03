function z = applyRegistration(dtrans,m,N)
%   z = applyRegistration(dtrans,m,z,N)

x=m.x;
y=m.y;
res = x(2)-x(1);

% to speed up/converve memory, we'll perform gridded interpolation over
% 2509x2509pixel subsets, giving 20x20 subsets = 400 total. Each subset
% overlaps by one pixel for merging.
sz = whos(m,'z'); sz = sz.size; % image dimensions info
int=1:2510:sz(2);
int=[int(1:end-1)',int(2:end)'];
int(1:end-1,2)=int(1:end-1,2)-1;

% need to add pixel buffer to int to account for shifts, scaling with dtrans
buff=max(ceil(abs(dtrans(2:3)./res)));
buff=[[0 buff];repmat([-buff buff],size(int,1)-2,1);[-buff 0]];

% initialize output
z = nan(sz,'single'); 

%nested interpolation loops
for i=1:length(int);
    
    % make subset row vector with buffer
    intRowBuff=int(i,1)+buff(i,1):int(i,2)+buff(i,2);
    
    for j=1:length(int);
        
        % make subset column vector with buffer
        intColBuff=int(j,1)+buff(j,1):int(j,2)+buff(j,2);
        
        % extract subset from unregistered tile with buffer and add offset
        z0 = m.z(intRowBuff,intColBuff) + dtrans(1);
        
        % apply coregistration cluster mask
        z0(~N(intRowBuff,intColBuff)) = NaN;
     
        if ~any(~isnan(z0(:))); continue; end
        
        % make unbuffered row and column vectors to insert into new
        % array
        intRow=int(i,1):int(i,2);
        intCol=int(j,1):int(j,2);
        
        % apply horizontal offsets and interpolate
        z(intRow,intCol)=interp2(x(intColBuff)+dtrans(2),...
                                 y(intRowBuff)+dtrans(3), z0,...
                                 x(intCol),y(intRow),'*linear');
        
        
    end
end