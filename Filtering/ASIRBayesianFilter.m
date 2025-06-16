function x = ASIRBayesianFilter(f, h, x0, P0, y, R, Q, Ns, Cx, Cy)
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
mu = zeros(N, Ns);
x(:, 1) = x0(:);
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
            PP = @(op, xt, xt_1) ft(op, xt, xt_1, Qt);
        else
            PP = @(op, xt, xt_1) ft(op, xt, xt_1, Qt, Cx_t);
        end
        mu(:, i) = PP('DS', 1, particles(:, i));
        if isempty(Cy_t)
            PL = @(yt, xt) ht(yt, xt, Rt);
        else
            PL = @(yt, xt) ht(yt, xt, Rt, Cy_t);
        end
        W(i) = PL(yt, mu(:, i)) * W(i);
    end
    W = W / sum(W);
    
    % Resampling
    index = randsample(1 : Ns, Ns, true, W);
    for i = 1 : Ns
        particles(:, i) = PP('DS', 1, particles(:, index(i)));
        W(i) = PL(yt,  particles(:, i)) / PL(yt,  mu(:, index(i)));
    end
    W =  W / sum(W);
    x(:, t) = sum(particles * diag(W), 2);
end