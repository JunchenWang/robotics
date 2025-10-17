function fX = AXXB_ndq(fA, fB)
fX0 = AXXB_de(fA, fB);
qd0 = Frame2qd(fX0);
f = @(qd) objFunc(qd, fA, fB);
qd = LevenbergMarquardt(f, qd0);
norm_qd(qd)
fX = qd2Frame(qd);

function [f J] = objFunc(qd, fA, fB)
n = size(fA, 2);
J = zeros(6 * n + 2, 8);
q = qd(1:4);
qe = qd(5:8);
lambda = 1e3;
S = zeros(6 * n, 8);
for i = 1 : n
    A = Frame_T_Converter(fA(:, i));
    B = Frame_T_Converter(fB(:, i));
    [thetaA dA lA mA] = ScrewAxis(A(1 : 3, 1 : 3), A(1 : 3, 4));
    [thetaB dB lB mB] = ScrewAxis(B(1 : 3, 1 : 3), B(1 : 3, 4));
    S(6 * i - 5 : 6 * i, :) = [lA - lB, SkewMatrix(lA + lB), zeros(3, 1), zeros(3)
                               mA - mB, SkewMatrix(mA + mB), lA - lB, SkewMatrix(lA + lB)];                      
end
J(1 : 6 * n, :) = S;
J(end - 1, :) = -2 * lambda * [q', 0 0 0 0];
J(end, :) = lambda * [qe', q'];
f =  [S * qd; lambda * (1 - dot(q, q)); lambda * dot(q, qe)];
