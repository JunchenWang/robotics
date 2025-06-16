function x = MH_MCMC(p, q, x0, N)
%ref: Understanding the Metropolis-Hastings Algorithm
%This is a general framework of Metropolis-Hastings algorithm
if nargin < 4
    N = 100;
end
x = repmat(x0, 1, N + 1);
for i = 2 : N + 1
    x_pre = x(:, i - 1);
    y = q('DS', 1, x_pre);
    u = rand;
    alpha = p(y) * q('v', x_pre, y) / (p(x_pre) * q('v', y, x_pre));
    if isnan(alpha)
        alpha = 1;
    else
        alpha = min(1, alpha);
    end
    if u <= alpha
        x(:, i) = y;
    else
        x(:, i) = x_pre;
    end
end

