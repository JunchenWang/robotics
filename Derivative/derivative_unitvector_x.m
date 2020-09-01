function J = derivative_unitvector_x(x)
J = -SkewMatrix(x)^2 / norm(x)^3;