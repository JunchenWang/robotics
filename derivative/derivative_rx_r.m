function J = derivative_rx_r(r, x)
r = r(:);
x = x(:);
theta = norm(r);
if (theta == 0)
    J = -SkewMatrix(x);
    return;
end
alpha = sin(theta) / theta;
beta = (1 - cos(theta)) / theta^2;
gama = (cos(theta) - alpha) / theta^2;
delta = (alpha - 2 * beta) / theta^2;
rrt = r * r';
SrSx = x * r' - dot(x, r) * eye(3);
J = -SkewMatrix(x) * (gama * rrt - beta * SkewMatrix(r) + alpha * eye(3))...
    - SrSx * (delta * rrt + 2 * beta * eye(3));
