function [Jf Jx] = derivative_fx(f, x)
f = f(:);
r = f(1 : 3);
x = reshape(x, 3, []);
N = size(x, 2);
Jf = zeros(3 * N, 6);
Jx = RotationByAxisAngleRep(r);
for i = 1 : N
    Jf(3 * i - 2 : 3 * i, :) = [derivative_rx_r(r, x(:, i)) eye(3)];
end