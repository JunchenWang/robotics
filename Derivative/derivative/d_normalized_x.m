function df = d_normalized_x(x, dx)
nx = norm(x);
df = (dx * nx - x * d_norm_x(x, dx)) / nx^2;