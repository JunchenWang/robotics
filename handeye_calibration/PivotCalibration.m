function [tip w] = PivotCalibration(r, t)
[tip w] = PivotCalibrationDLT(r, t);
x0 = [r(:); tip; w];
f = @(x) objFun(x, r, t);
x = LevenbergMarquardt(f, x0, [], [], [], []);
tip = x(end - 5 : end - 3);
w = x(end - 2 : end);

function [f J] = objFun(x, r, t)
N = size(r, 2);
f = zeros(6 * N, 1);
J = zeros(6 * N, 3 * N + 6);
tip = x(end - 5 : end - 3);
w = x(end - 2 : end);
for i = 1 : N
    ri = x(3 * i - 2 : 3 * i);
    
    f(6 * i - 5 : 6 * i) = [ri - r(:, i); w - ApplyRotationToPoint(ri, tip) - t(:, i)];
    
    J(6 * i - 5 : 6 * i, [3 * i - 2 : 3 * i, end - 5 : end]) ...
        = [eye(3), zeros(3, 6); ...
        -derivative_rx_r(ri, tip), -RotationByAxisAngleRep(ri), eye(3)];
end