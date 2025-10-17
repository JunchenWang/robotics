function [tip w] = PivotCalibrationDLT(r, t)
N = size(r, 2);
A = zeros(3 * N, 6);
b = zeros(3 * N, 1);
for i = 1 : N
    ri = r(:, i);
    ti = t(:, i);
    Ri = RotationByAxisAngleRep(ri);
    A(3 * i - 2 : 3 * i, :) = [Ri, -eye(3)];
    b(3 * i - 2 : 3 * i) = -ti;
end
x = A \ b;
tip = x(1:3);
w = x(4 : end);



