function fX = NonlinearAXXB(fA, fB, lambda1, lambda2, lambda)
if nargin < 3
    lambda1 = 1;
end
if nargin < 4
    lambda2 = 1;
end
if nargin < 5
    lambda = sqrt(2) * 1e3;
end
f = @(x) objFunc(x(1 : 4), x(5 : 7), fA, fB, lambda1, lambda2, lambda);
fX0 = AXXB_de(fA, fB);
q0 = r_q_converter(fX0(1:3));
t0 = fX0(4 : 6);
x0 = [q0; t0];
% opt = optimset(@lsqnonlin);
% opt = optimset(opt, 'Display',  'iter', 'Algorithm', 'levenberg-marquardt','Jacobian', 'on', 'TolFun', 1e-12, 'TolX', 1e-12, 'DerivativeCheck', 'on');
% [x,err,RESIDUAL,EXITFLAG,OUTPUT] = lsqnonlin(f, x0, [], [], opt);
x = LevenbergMarquardt(f, x0, 1e-12, 1e-12, 1e-12, [], 0);
fX = [r_q_converter(x(1 : 4)); x(5 : 7)];

function [f J] = objFunc(q, t, fA, fB, lambda1, lambda2, lambda)
N = size(fA, 2);
f = zeros(8 * N + 1, 1);
J = zeros(8 * N + 1, 7);

for i = 1 : N
    fAi = fA(:, i);
    fBi = fB(:, i);
    tA = fAi(4 : 6);
    nA = ExtractKTheta(fAi);
    RA = RotationByAxisAngleRep(fAi(1 : 3));
    tB = fBi(4 : 6);
    nB = ExtractKTheta(fBi);
    tem1 = ConcatenateQuaternions([0;nA], q);
    tem2 = ConcatenateQuaternions(q, [0;nB]);
    f1 = tem1 - tem2;
    
    tem1 = ConcatenateQuaternions(q, [0;tB]);
    tem2 = ConcatenateQuaternions([0; (RA - eye(3)) * t], q);
    tem3 = ConcatenateQuaternions([0; tA], q);
    f2 = tem1 - tem2 - tem3;
    f(8 * i - 7 : 8 * i) = [lambda1 * f1; lambda2 * f2];
    
    [J1 J2] = derivative_q2starq1([0;nA], q);
    [J3 J4] = derivative_q2starq1(q, [0;nB]);
    Jq = J2 - J3;
    J(8 * i - 7 : 8 * i - 4, :) = lambda1 * [Jq, zeros(4, 3)]; 
    [J1 J2] = derivative_q2starq1(q, [0;tB]);
    [J3 J4] = derivative_q2starq1([0; (RA - eye(3)) * t], q);
    [J5 J6] = derivative_q2starq1([0;tA], q);
    Jq = J1 - J4 - J6;
    Jt = -J3 * [zeros(1, 3); RA - eye(3)];
    J(8 * i - 3 : 8 * i, :) = lambda2 * [Jq, Jt]; 
end
f(end) = lambda * (1 - dot(q, q));
J(end, :) = [-2 * lambda * q', 0 0 0];