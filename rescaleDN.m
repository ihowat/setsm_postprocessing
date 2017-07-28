function dn=rescaleDN(or,dnmax)
% RESCALEDN rescale digitial numbers to new maximum
%
% dn=rescaleDN(dn,dnmax)rescales the digital number image dn to the new
% maximum in dnmax.
%
% Ian Howat, ihowat@gmail.com
% 24-Jul-2017 15:50:25

% Set the minimum and maximum values of this scale. We use a fixed scale
% because this is what all data is scaled to after application of
% wv_correct regardless of actual min or max.
ormin = 0; 
ormax = 32767;

% Set the new minimum and maximum. The dnmin is zero because nodata is
% apparently used in the scaling.
dnmin = 0;% single(stats.min(stats_ind));
dnmax = single(dnmax);

% rescale back to original dn
dn = dnmin + (dnmax-dnmin).*(single(or) - ormin)./(ormax-ormin);

