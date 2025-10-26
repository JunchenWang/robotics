function J = derivative_Nx(x)
J = x(:)' / norm(x);