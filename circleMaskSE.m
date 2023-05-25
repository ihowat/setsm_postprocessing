function SE = circleMaskSE(radius)

n_approx = 4;

SE = strel("disk", radius, n_approx);
