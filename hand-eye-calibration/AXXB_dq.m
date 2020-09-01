function x = AXXB_dq(fA, fB)
n = size(fA, 2);
S = zeros(6 * n, 8);
for i = 1 : n
    A = Frame_T_Converter(fA(:, i));
    B = Frame_T_Converter(fB(:, i));
    [thetaA dA lA mA] = ScrewAxis(A(1 : 3, 1 : 3), A(1 : 3, 4));
    [thetaB dB lB mB] = ScrewAxis(B(1 : 3, 1 : 3), B(1 : 3, 4));
    S(6 * i - 5 : 6 * i, :) = [lA - lB, SkewMatrix(lA + lB), zeros(3, 1), zeros(3)
                               mA - mB, SkewMatrix(mA + mB), lA - lB, SkewMatrix(lA + lB)];
                           
end
[U D V] = svd(S, 0);
u1 = V(1 : 4, 7);
v1 = V(5 : 8, 7);
u2 = V(1 : 4, 8);
v2 = V(5 : 8, 8);

u11 = dot(u1, u1);
u12 = dot(u1, u2);
u22 = dot(u2, u2);
u1v1 = dot(u1, v1);
u1v2 = dot(u1, v2);
u2v1 = dot(u2, v1);
u2v2 = dot(u2, v2);

if abs(u1v1) > abs(u2v2)
    a = u1v1;
    b = u1v2 + u2v1;
    c = u2v2;
    s1 = (-b + sqrt(b^2 - 4 * a *c)) / (2 * a);
    s2 = (-b - sqrt(b^2 - 4 * a *c)) / (2 * a);
    test1 = s1^2 * u11 + 2 * s1 * u12 + u22;
    test2 = s2^2 * u11 + 2 * s2 * u12 + u22;
    if  test1 > test2
        s = s1;
        test = test1;
    else
        s = s2;
        test = test2;
    end
    lambda2 = sqrt(1 / test);
    lambda1 = s * lambda2;
else
    c = u1v1;
    b = u1v2 + u2v1;
    a = u2v2;
    s1 = (-b + sqrt(b^2 - 4 * a *c)) / (2 * a);
    s2 = (-b - sqrt(b^2 - 4 * a *c)) / (2 * a);
    test1 = s1^2 * u22 + 2 * s1 * u12 + u11;
    test2 = s2^2 * u22 + 2 * s2 * u12 + u11;
    if  test1 > test2
        s = s1;
        test = test1;
    else
        s = s2;
        test = test2;
    end
    lambda1 = sqrt(1 / test);
    lambda2 = s * lambda1;
end
qdx = lambda1 * [u1; v1] + lambda2 * [u2; v2];
[R t] = transform_qd(qdx);
x = Frame_T_Converter([R t; 0 0 0 1]);