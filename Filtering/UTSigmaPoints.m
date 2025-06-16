function [Px W WP] = UTSigmaPoints(mu, sigma, alpha, beta, k)
%Scale Unscented Transform
mu = mu(:);
N = length(mu);
if nargin < 3
    alpha = 1;
end
if nargin < 4
    beta = 0;
end
if nargin < 5
    k = 3 - N;
end
Px = zeros(N, 2 * N + 1);

w0 = k / (N + k);
W = [1 / (2 * (N + k)) * ones(2 * N, 1); w0] / alpha^2;
W(end) = W(end) + (1 - 1 / alpha^2);
WP = W;
WP(end) = WP(end) + 1 - alpha^2 + beta;
[R p] = chol((N + k) * sigma);
if p ~= 0
    error('chol error!');
end
for i = 1 : N
    Px(:, i) = mu + alpha * R(i, :)';
    Px(:, N + i) = mu - alpha * R(i, :)';
end
Px(:, 2 * N + 1) = mu;
