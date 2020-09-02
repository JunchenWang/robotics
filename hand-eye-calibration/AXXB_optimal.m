function fX = AXXB_optimal(fA, fB)
% fX0 = NonlinearAXXB(fA, fB);
fX0 = AXXB_de(fA, fB);
x0 = [fB(:); fX0];
f = @(x) Func(x, fA, fB);
opt = optimset(@lsqnonlin);
opt = optimset(opt, 'Display',  'off', 'Algorithm', 'levenberg-marquardt','Jacobian', 'on', 'TolFun', 1e-12, 'TolX', 1e-12, 'DerivativeCheck', 'off');
[x,err,RESIDUAL,EXITFLAG,OUTPUT] = lsqnonlin(f, x0, [], [], opt);
% x = LevenbergMarquardt(f, x0, 1e-12, 1e-12, 1e-12, [], 0);
fX = x(end - 5 : end);


function [f, J] = Func(x, fA, fB)
n = size(fA, 2);
f = zeros(12 * n, 1);
fX = x(end - 5 : end);
J = zeros(12 * n, 6 * n + 6);
lambda = 15;
weight = diag([lambda lambda lambda 1 1 1]);
for i = 1 : n
    fB_bar = x(6 * i - 5 : 6 * i);
    fX_ = InverseFrame(fX);
    fxbx_ = ConcatenateFrame(fX, ConcatenateFrame(fB_bar, fX_));
    f(12 * i - 11 : 12 * i - 6) = weight * (fB_bar - fB(:, i));
    J(12 * i - 11 : 12 * i - 6, 6 * i - 5 : 6 * i) = weight;
    f(12 * i - 5 : 12 * i) = weight * (fxbx_ - fA(:, i));
    [Jfx, Jfb] = derivative_fxfbfx_(fX, fB_bar);
    J(12 * i - 5 : 12 * i, [6 * i - 5 : 6 * i, end - 5 : end]) = weight * [Jfb, Jfx];
end
