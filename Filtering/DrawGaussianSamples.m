function x = DrawGaussianSamples(f, mu, sigma, N)
M = length(mu);
mu = mu(:);
x = zeros(M, N);
s = repmat(mu, 1, N) + chol(sigma) * randn(M, N);
for i = 1 : N
    x(:, i) = f(s(:, i));
end