function [Y, YTr] = regressor_matrix(robot, q, qd, a, v, r)
% paper: On the closed form computation of the dynamic matrices and their differentiations
n = robot.dof;
A = robot.A;
M = robot.M;
J = zeros(6, n, n);
dJ = zeros(6, n, n);
gravity = [0;0;0;-robot.gravity'];
YTr = zeros(10 * n, 1);
Y = zeros(n, 10 * n);
for i = 1 : n
    T = eye(4);
    dT = zeros(4,4);% T对时间的导数a
    for j = i : -1: 1
        [dAdT, AdT] = derivative_adjoint_T(T, dT);
        J(:,j,i) = AdT * A(j,:)';
        dJ(:,j,i)= dAdT * A(j,:)';
        tform = exp_twist(-A(j,:) * q(j)) * tform_inv(M(:,:,j));
        dT = (dT  + T * se_twist(-A(j,:)) * qd(j)) * tform;
        T = T * tform;
    end
    Jk = J(:,:,i);
    dJk = dJ(:,:,i);
    Adk = adjoint_T(T);
    Vk = Jk * qd;
    adk = adjoint_twist(Vk);
    rk = Adk * gravity;
    alpha = Jk * a + adk * Jk * v + dJk * v +  rk;
    Y(:, (i - 1) * 10 + 1 : 10 * i) = Jk' * (A_matrix(alpha) - adk' * A_matrix(Jk * v));
    YTr((i - 1) * 10 + 1 : 10 * i) = Y(:, (i - 1) * 10 + 1 : 10 * i)' * r;
end
end

function A = A_matrix(V)
w = V(1:3);
v = V(4:6);
Aw = [zeros(3,1), -so_w(v), diag(w), [w(2), w(3), 0; w(1), 0, w(3); 0, w(1), w(2)]];
Av = [v, so_w(w), zeros(3,6)];
A = [Aw; Av];
end

