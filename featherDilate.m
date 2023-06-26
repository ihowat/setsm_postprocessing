function A = featherDilate(I,npix)

if npix <= 1
    A = I;
    return;
end

if mod(npix, 2) == 0
    npix = npix + 1;
end
npix_half = (npix - 1) / 2;

if npix_half == 1
    npix_half = 0;
end

circle_buff_mask = circleMaskSE(npix_half+1);
circle_feather_mask = circleFeatherMask(npix+2, 0);
circle_feather_mask_ones = (circle_feather_mask ~= 0);

circle_feather_mask = 1 - circle_feather_mask;
circle_feather_mask(circle_feather_mask == 1) = 0;
circle_feather_mask(circle_feather_mask == 0 & circle_feather_mask_ones == 1) = 1;

I = padarray(I, [npix_half, npix_half], 0);
A = imdilate(I, circle_buff_mask);

% mask = conv2(double(mask), circle_feather_mask, 'same') ./ conv2(ones(size(mask),'double'), circle_feather_mask, 'same');

% mask = conv2(double(mask), circle_feather_mask, 'same') ./ conv2(ones(size(mask),'double'), circle_feather_mask_ones, 'same');
% mask = mask / (sum(circle_feather_mask, 'all') / sum(circle_feather_mask_ones, 'all'));

A = conv2(A, circle_feather_mask_ones, 'same');
A = A / sum(circle_feather_mask_ones, 'all');

A = A(npix_half+1:end-npix_half, npix_half+1:end-npix_half);
