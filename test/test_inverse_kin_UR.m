function   [err, failure, cnt, t] = test_inverse_kin_UR(N)
if nargin < 1
    N = 10000;
end
% a2 a3 在dh表里为负，这里给长度即可
param = [0.425, 0.3922, 0.1625, 0.1333, 0.0997, 0.0996];
tol = [1e-5, 1e-5];
% cfg = [1,1,1];
n = 6;
err = zeros(N, 2);
failure = zeros(1000, n);
cnt = 0;
tic;
for i =1 : N
    q = rand(n, 1) * 2 * pi - pi;
    q(5) = 0;
    Td = forward_kin_UR(q, param);
    [angles, flag] = inverse_kin_UR_near(Td, param, q, tol);
   
    T = forward_kin_UR(angles, param);
    err(i, 1) = norm(logR(T(1:3,1:3)'*Td(1:3,1:3)));
    err(i, 2) = norm(T(1:3,4)- Td(1:3,4));
    if ~flag
        cnt = cnt + 1;
        failure(cnt, :) = q;
    end
end
failure(cnt + 1:end, :) = [];
t = toc;
end

