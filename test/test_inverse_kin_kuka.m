function   [err, failure, cnt, t] = test_inverse_kin_kuka(N)
if nargin < 1
    N = 10000;
end
tol = [1e-5, 1e-5];
kesai = 0;
n = 7;
err = zeros(N, 3);
failure = zeros(1000, n);
cnt = 0;
limit = pi * 0.8;
tic;
for i =1 : N
    q = rand(n, 1) * 2 * limit - limit;
    q(3) = 0;
    Td = forward_kin_kuka(q);
    kesai = cal_kuka_kesai(q);
    [angles, flag] = inverse_kin_kuka_kesai_near(Td, kesai, q, tol);
   
    T = forward_kin_kuka(angles);
    err(i, 1) = norm(logR(T(1:3,1:3)'*Td(1:3,1:3)));
    err(i, 2) = norm(T(1:3,4)- Td(1:3,4));
    err(i, 3) = norm(angles - q);
    % if err(i,3) > 1e-2
    %     cnt = cnt + 1;
    %     failure(cnt, :) = q;
    % end
    if ~flag
        cnt = cnt + 1;
        failure(cnt, :) = q;
    end
end
failure(cnt + 1:end, :) = [];
t = toc;
disp(mean(err, 1));
end


