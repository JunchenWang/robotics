function df = d_norm_x(x, dx)
df =  x'* dx / sqrt(x'*x);