function [x P] = EKalmanFilter(f, h, x0, P0, y, R, Q, Cx, Cy)
% display('EK');
error(nargchk(7, 9, nargin));
if nargin < 8
    Cx = [];
end
if nargin < 9
    Cy = [];
end
KFPredictionFunc = @(ft, ht, Rt, Qt, mu_x, sigma_x, Cx_t, Cy_t) ...
    KalmanPrediction(ft, ht, Rt, Qt, mu_x, sigma_x, Cx_t, Cy_t, 'E');
[x P] = GeneralKalmanFilter(f, h, x0, P0, y, R, Q, KFPredictionFunc, Cx, Cy);
