function fX = AXXB_dq2(fA, fB)
n = size(fA, 2);
K1 = zeros(3 * n, 4);
K2 = zeros(3 * n, 4);
for i = 1 : n
    A = Frame_T_Converter(fA(:, i));
    B = Frame_T_Converter(fB(:, i));
    [thetaA dA lA mA] = ScrewAxis(A(1 : 3, 1 : 3), A(1 : 3, 4));
    [thetaB dB lB mB] = ScrewAxis(B(1 : 3, 1 : 3), B(1 : 3, 4));
    K1(3 * i - 2 : 3 * i, :) = [lA - lB, SkewMatrix(lA + lB)]; 
    K2(3 * i - 2 : 3 * i, :) = [mA - mB, SkewMatrix(mA + mB)];   
end
[U D V] = svd(K1, 0);
q = V(:, end);
basis = null(q');
b = K2 * q;
lambda = K1 * basis \ -b;
qe = basis * lambda;
qd = [q; qe];
fX = qd2Frame(qd);


