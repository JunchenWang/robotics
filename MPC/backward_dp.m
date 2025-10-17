function [K, U] = backward_dp(A, B, Q, R, Pf, N, x0)
% x[n + 1] = A * x[n] + B * u[n]
% minimize sum 1/2 * (x[k]'Qx[k] + u[k]'Ru[k]) + 1/2 * x[n]'Qfx[n]
n = size(A, 2);
m = size(B, 2);
K = zeros(m, n, N);
PI = Pf;
U = zeros(m, N);
for k = N : -1 : 1
    K(:,:,k) = (B' * PI * B + R) \ (B' * PI * A);
    PI = Q + A' * PI * A - A' * PI * B * K(:,:,k);
end
x = x0;
for k = 1 : N
    U(:, k) = -K(:,:,k) * x;
    x = (A - B * K(:,:,k)) * x;
end
