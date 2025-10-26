function [G SigmaG] = MonteCarlo(n, g, p)
if nargin == 1
    g = @(x) x;
    p = [];
elseif nargin == 2
    p = [];
end
if isempty(p)
    x = randn(1, n);
else
    x = p('DS', n);
end
v = g(x(:, 1));
mu = v;
ss = v * v';
for i = 2 : n
    v = g(x(:, i));
    mu = 1 / i * ((i - 1) * mu + v);
    ss = ss + v * v';
%     sigma = 1 / i * (1 / (i - 1) * ss - i / (i - 1) * (mu * mu'));
end
sigma = 1 / n * (1 / (n - 1) * ss - n / (n - 1) * (mu * mu'));
G = mu;
SigmaG = sigma;