function [angles, flag] = inverse_kin_general(robot, Td, ref, tol, max_iter)
if nargin < 4
    tol = [1e-4, 1e-4];
end
if nargin < 5
    max_iter = 100;
end
flag = 0;
alpha = 1;
lambda = 0.01;
max_step = 0.1;
q = ref(:);
Rd = Td(1:3,1:3);
pd = Td(1:3, 4);
for i = 1 : max_iter
    [J, T] = jacobian_matrix(robot, q);
    R = T(1:3, 1:3);
    p = T(1:3, 4);
    xe = [logR(R' * Rd)'; R' * (pd - p)];
    if norm(xe(1:3)) < tol(1) && norm(xe(4:6)) < tol(2)
        angles = q;
        flag = 1;
        return;
    end
    qe = damping_least_square(J, xe, lambda, 0);
    qe_norm = norm(qe);
    if qe_norm > max_step
        qe = qe / qe_norm * max_step;
    end
    q = q + alpha * qe;
end
angles = q;
