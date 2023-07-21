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

buff_mask = circleMaskSE(npix_half+1);
feather_mask = circleMaskDiameter(npix+1);
% buff_mask = ones(npix);
% feather_mask = ones(npix);

I = padarray(I, [npix_half, npix_half], 0);
A = imdilate(I, buff_mask);

A = conv2(A, feather_mask, 'same');
A = A / sum(feather_mask, 'all');

A = A(npix_half+1:end-npix_half, npix_half+1:end-npix_half);
