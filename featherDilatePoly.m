function A = featherDilatePoly(xref,yref,I,xi,yi,buff_d)

ps = polyshape(xi, yi);
ps_buff = polybuffer(ps, buff_d);
%ps_buff = polybuffer(ps, buff_d, 'JointType','miter');

dx = abs(xref(2)-xref(1));
buff_px = ceil(buff_d / dx);

BW0 =  roipoly(xref,yref,I,     ps.Vertices(:,1),     ps.Vertices(:,2));
BW1 = ~roipoly(xref,yref,I,ps_buff.Vertices(:,1),ps_buff.Vertices(:,2));

A0 = featherDilate(BW0, buff_px);
A1 = 1 - featherDilate(BW1, buff_px);

A = (A0 + A1) / 2;

% Fix edge case issues
A(BW0) = 1;
A(BW1) = 0;
