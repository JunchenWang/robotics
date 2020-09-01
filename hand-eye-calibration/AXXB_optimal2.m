function fX = AXXB_optimal2(fA, fB)
qd_x0 = Frame2qd(AXXB_de(fA, fB));
qd_b = Frame2qd(fB);
x0 = [qd_b(:); qd_x0];
f = @(x) objFunc(x, fA, fB);
opt = optimset(@lsqnonlin);
opt = optimset(opt, 'Display',  'iter', 'Algorithm', 'levenberg-marquardt','Jacobian', 'off', 'TolFun', 1e-12, 'TolX', 1e-12, 'DerivativeCheck', 'off');
[x,err,RESIDUAL,EXITFLAG,OUTPUT] = lsqnonlin(f, x0, [], [], opt);
fX = qd2Frame(x(end - 7 : end));

function [f J] = objFunc(x, fA, fB)
n = size(fA, 2);
f = zeros(16 * n + 2, 1);
J = zeros(16 * n, 8 * n + 8);
qd_x = x(end - 7 : end);
lambda = 1e4;
for i = 1 : n
    qd_b = x(8 * i - 7 : 8 * i);
    qd_x_ = Conjugate_qd(qd_x);
    fxbx_ = ConcatenateDualQuaternions(qd_x, qd_b, qd_x_);
    f(16 * i - 15 : 16 * i - 8) = qd_b - Frame2qd(fB(:, i));
    qd_a = Frame2qd(fA(:, i));
    if fxbx_(1) < 0;
        qd_a = -qd_a;
    end
    f(16 * i - 7 : 16 * i) = fxbx_ - qd_a;
%     J(12 * i - 11 : 12 * i - 6, 6 * i - 5 : 6 * i) = eye(6);
%     [J1 J2] = derivative_f2f1(fX, ConcatenateFrame(fB_bar, fX_));
%     [J3 J4] = derivative_f2f1(fB_bar, fX_);
%     J5 = derivative_finv(fX);
%     J(12 * i - 5 : 12 * i, [6 * i - 1 : 6 * i, end - 5 :end]) = [J2 * J3, ];
end
 f(end - 1 : end) = [lambda * (1 - dot(qd_x(1:4), qd_x(1 : 4)))
                     lambda * dot(qd_x(1:4), qd_x(5 : 8))];