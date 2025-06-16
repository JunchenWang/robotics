function [mu_x sigma_x mu_y sigma_y sigma_xy] = UT_Cross(mu, sigma, f, h, alpha, beta, k)

error(nargchk(4, 7, nargin));
if nargin < 5
    alpha = 0.25; %1e-3;
end
if nargin < 6
    beta = 3; %2;
end
if nargin < 7
    k = 0;% 3 - length(mu);
end
[Px W WP] = UTSigmaPoints(mu, sigma, alpha, beta, k);

N = length(W);

zf0 = f(Px(:, 1));
Df = length(zf0);

zh0 = h(Px(:, 1));
Dh = length(zh0);

zf = zeros(Df, N);
zh = zeros(Dh, N);
zf(:, 1) = zf0;
zh(:, 1) = zh0;

sigma_x = zeros(Df);
sigma_y = zeros(Dh);
sigma_xy = zeros(Df, Dh);
mu_x = W(1) * zf(:, 1);
mu_y = W(1) * zh(:, 1);

for i = 2 : N
    zf(:, i) = f(Px(:, i));
    mu_x = mu_x +  W(i) * zf(:, i);
    
    zh(:, i) = h(Px(:, i));
    mu_y = mu_y +  W(i) * zh(:, i);
end

for i = 1 : N
    sigma_x = sigma_x +  WP(i) *  (zf(:, i) - mu_x) * (zf(:, i) - mu_x)';
    sigma_y = sigma_y +  WP(i) *  (zh(:, i) - mu_y) * (zh(:, i) - mu_y)';
    sigma_xy = sigma_xy +  WP(i) *  (zf(:, i) - mu_x) * (zh(:, i) - mu_y)';
end

