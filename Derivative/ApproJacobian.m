function J = ApproJacobian(f, x, step)
y = f(x);
M = length(y);
N = length(x);
J = zeros(M, N);
if nargin < 3
    step = 1e-4;
end

for j = 1 : N
    delta = zeros(N, 1);
    delta(j) = step;
    J(:, j) = (f(x + delta) - y) / step;
end

