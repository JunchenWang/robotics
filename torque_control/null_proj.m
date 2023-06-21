function Nx = null_proj(J, M, x)
Nx = x - pinv_J_x(J, M, J * x);