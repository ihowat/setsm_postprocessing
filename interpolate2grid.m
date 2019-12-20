function [z,mt,or] = interpolate2grid(x,y,z,mt,or,xi,yi)
% interpolate2grid interpolate the mask and ortho images to shifted grid
%
% [mi,oi] = interpolate2grid(m,o,x,y,dtrans)

% interpolate elevation to grid
z = interp2(x,y,z,xi(:)',yi(:),'*linear');

% interpolate the mask to the same grid
mt = interp2(x,y,single(mt),xi(:)',yi(:),'*nearest');

mt(isnan(mt)) = 0; % convert back to uint8
mt = logical(mt);


% interpolate ortho to same grid
or = single(or);
or(or==0) = NaN; % set border to NaN so wont be interpolated
try

 or = interp2(x,y,or,xi(:)',yi(:),'*cubic');

catch
  save interpolate2grid_dump.mat
 
  error('hit interpolate2grid error, dumped content to  interpolate2grid_dump.mat')

end
or(isnan(or)) = 0; % convert back to int16
or = int16(or);
