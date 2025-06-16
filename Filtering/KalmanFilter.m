function [x P] = KalmanFilter(f, h, x0, P0, y, R, Q, Cx, Cy)

error(nargchk(7, 9, nargin));
if nargin < 8
    Cx = [];
end
if nargin < 9
    Cy = [];
end
KFPredictionFunc = @(ft, ht, Rt, Qt, mu_x, sigma_x, Cx_t, Cy_t) ...
    KalmanPrediction(ft, ht, Rt, Qt, mu_x, sigma_x, Cx_t, Cy_t, 'L');
[x P] = GeneralKalmanFilter(f, h, x0, P0, y, R, Q, KFPredictionFunc, Cx, Cy);