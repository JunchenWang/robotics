function x = SIRBayesianFilter(f, h, x0, P0, y, R, Q, Ns, Cx, Cy)
if nargin < 8
    Ns = 50;
end
if nargin < 9
     Cx = [];
end
if nargin < 10
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
    index = randsample(1 : Ns, Ns, true, W);
    W = ones(1, Ns) / Ns;
    particles = particles(:, index);
    x(:, t) = sum(particles * diag(W), 2);
end