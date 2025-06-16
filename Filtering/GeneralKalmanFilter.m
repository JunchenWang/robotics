function [x P] = GeneralKalmanFilter(f, h, x0, P0, y, R, Q, KFPredictionFunc, Cx, Cy)

if nargin < 9
    Cx = [];
end
if nargin < 10
    Cy = [];
end

[T, f, h, y, R, Q, Cx, Cy, N] = FilterInputPreprocess(f, h, y, R, Q, Cx, Cy, x0);
x = zeros(N, T + 1);
P = zeros(N, N, T + 1);
x(:, 1) = x0;
P(:, :, 1) = P0;
for t = 2 : T + 1
    Cx_t = Cx{t - 1};
    Cy_t = Cy{t - 1};
    Rt = R{t - 1};
    Qt = Q{t - 1};
    ft = f{t - 1};
    ht = h{t - 1};
    if isa(Qt, 'function_handle')
        Qt = Qt(x(:, t - 1));
    end
    DimQ = size(Qt, 1);
    if isa(Rt, 'function_handle')
        if isempty(Cx_t)
            xt = ft(x(:, t - 1), zeros(DimQ, 1));
        else
            xt = ft(x(:, t - 1), zeros(DimQ, 1), Cx_t);
        end
        Rt = Rt(xt);
    end
    yt = y{t - 1};
    [mu_x_t sigma_x_t mu_y_t sigma_y_t sigma_xy_t] = ...
        KFPredictionFunc(ft, ht, Rt, Qt, x(:, t - 1), P(:, :, t - 1), Cx_t, Cy_t);
    [x(:, t) P(:, :, t)] = KalmanUpdation(mu_x_t, sigma_x_t, mu_y_t, sigma_y_t, sigma_xy_t, yt);
end