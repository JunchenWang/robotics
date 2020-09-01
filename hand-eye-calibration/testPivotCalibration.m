function [err1 err2] = testPivotCalibration
N = 100;
r = zeros(3, N);
t = r;
tip = [10, 25, -160]';
w = [20, 40, -500]';
for i = 1 : N
ri = rand(3, 1);
ri = ri / norm(ri);
thetai = rand * pi;
ri = thetai * ri;
ti = w - ApplyRotationToPoint(ri, tip);
t(:, i) = ti;
r(:, i) = ri;
end

M = 1000;
err1 = zeros(6, M);
err2 = err1;
for i = 1 : M
    r_n = r + 0.05 * randn(3, N);
    t_n = t + 0.5 * randn(3, N);
    
    [tip_e w_e] = PivotCalibrationDLT(r_n, t_n);
    err1(:,i) = [tip_e; w_e] - [tip; w];
    [tip_e w_e] = PivotCalibration(r_n, t_n);
    err2(:,i) = [tip_e; w_e] - [tip; w];
end
