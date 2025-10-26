function x = RegularizedBayesianFilter(f, h, x0, P0, y, R, Q, Ns, Cx, Cy)
if nargin == 7
    Ns = 50;
    Cx = [];
    Cy = [];
end
if nargin == 8
     Cx = [];
     Cy = [];
end
if nargin == 9
     Cy = [];
end
[T, f, h, y, R, Q, Cx, Cy, N] = FilterInputPreprocess(f, h, y, R, Q, Cx, Cy, x0);
% P = zeros(N, Ns, T + 1);
% P(:, :, 1) = P0('DS', Ns);
particles = P0('DS', Ns);
W = ones(1, Ns) / Ns;
x = zeros(N, T + 1);
x(:, 1) = x0;
for t = 2 : T + 1
    Cx_t = Cx{t - 1};
    Cy_t = Cy{t - 1};
    Rt = R{t - 1};
    Qt = Q{t - 1};
    ft = f{t - 1};
    ht = h{t - 1};
    yt = y{t - 1};
    for i = 1 : Ns
        if isempty(Cx_t)
            PP = @(op, xt) ft(op, xt, particles(:, i), Qt);
        else
            PP = @(op, xt) ft(op, xt, particles(:, i), Qt, Cx_t);
        end
        particles(:, i) = PP('DS', 1);
        if isempty(Cy_t)
            PL = @(yt) ht(yt, particles(:, i), Rt);
        else
            PL = @(yt) ht(yt, particles(:, i), Rt, Cy_t);
        end
         W(i) = PL(yt);
    end
    W = W / sum(W);
    % Resampling
    mu = sum(particles * diag(W), 2);
    sigma = zeros(N);
    for i = 1 : Ns
        sigma = sigma + W(i) * (particles(:, i) - mu) * (particles(:, i) - mu)';
    end
    D = chol(sigma);
    D = D';
    index = randsample(1 : Ns, Ns, true, W);
    W = ones(1, Ns) / Ns;
    particles = particles(:, index);
    A = (8 / Hypersphere(N) * (N + 4) * (2 * sqrt(pi))^N)^(1 / (N + 4));
    hopt = A * Ns^(1 / (N + 4));
    
    for i = 1 : Ns
%         s = DrawSamples(@(x) kernel(x, N), 1, -1, 1);
        s = DrawSamples;
        particles(:, i) = particles(:, i) + hopt * D * s;
        
    end
    x(:, t) = sum(particles * diag(W), 2);
end

function p = kernel(x, nx)
N = size(x, 2);
p = zeros(N, 1);
for i = 1 : N
    if norm(x(:, i)) <= 1
        p(i) = (nx + 2) / (2 * Hypersphere(nx)) * (1 - norm(x(:, i))^2);
    else
        p(i) = 0;
    end
end

function x = DrawSamples
p = rand;
x = roots([-0.25 0 0.75 0.5 - p]);
x = x(3);

