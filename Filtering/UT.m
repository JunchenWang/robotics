function [mu_y sigma_y] = UT(mu_x, sigma_x, f, alpha, beta, k)

error(nargchk(3, 6, nargin));
if nargin < 4
    alpha = 0.25;%1e-3;
end
if nargin < 5
    beta = 3; %2;
end
if nargin < 6
    k = 0; %3 - length(mu_x);
end
[Px W WP] = UTSigmaPoints(mu_x, sigma_x, alpha, beta, k);
N = length(W);
y0 = f(Px(:, 1));
M = length(y0);
y = zeros(M, N);
y(:, 1) = y0;
sigma_y = zeros(M);
mu_y = W(1) * y(:, 1);
for i = 2 : N
    y(:, i) = f(Px(:, i));
    mu_y = mu_y +  W(i) * y(:, i);
end
for i = 1 : N
    sigma_y = sigma_y +  WP(i) *  (y(:, i) - mu_y) * (y(:, i) - mu_y)';
end

%%
% n = (length(W) - 1) / 2; 
% if nargin == 3
%     alpha = 1 / n;
%     mu = alpha^2;
% end
% if nargin == 4
%     mu = alpha^2;
% end
% g = @(x) (f(mu_x + alpha * (x - mu_x)) - f(mu_x)) / mu + f(mu_x);
% N = size(Px, 2);
% y0 = g(Px(:, 1));
% M = length(y0);
% y = zeros(M, N);
% y(:, 1) = y0;
% sigma_y = zeros(M);
% mu_y = W(1) * y(:, 1);
% for i = 2 : N
%     y(:, i) = g(Px(:, i));
%     mu_y = mu_y +  W(i) * y(:, i);
% end
% for i = 1 : N
%     sigma_y = sigma_y +  W(i) *  (y(:, i) - mu_y) * (y(:, i) - mu_y)';
% end
% sigma_y = mu * sigma_y;