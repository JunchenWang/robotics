function [Ym, Yg] = regressor_for_mg(robot, q)
% Ym: n(n+1)/2 by 10n, Ym theta yields low triangle of Mass. Yg: n by 10n
% theta: 
% paper: Dynamic Identification of the Franka Emika Panda Robot With Retrieval of Feasible Parameters Using Penalty-Based Optimization
n = robot.dof;
Ym2 = zeros(n, 10 * n, n);
Ym = zeros(n * (n + 1) / 2, 10 * n);
Yg = zeros(n, 10 * n);
A = robot.A;
M = robot.M;
J = zeros(6, n, n);
T_inv = zeros(4,4,n);
gravity = [0;0;0;-robot.gravity'];
for i = n : -1 : 1
   T = eye(4);
   for j = i : -1: 1
       J(:,j,i) = adjoint_T(T) * A(j,:)';
       T = T * exp_twist(-A(j,:) * q(j)) * tform_inv(M(:,:,j));
   end
   T_inv(:,:,i) = T;
end

for i = 1 : n
    Jk = J(:,:,i);
    rk = adjoint_T(T_inv(:,:,i)) * gravity;
    Yg(:, (i - 1) * 10 + 1 : 10 * i) = Jk' * A_matrix(rk);
    for j = 1 : n
        Ym2(:, (i - 1) * 10 + 1 : 10 * i, j) = Jk' * A_matrix(Jk(:, j) +  rk) - Yg(:, (i - 1) * 10 + 1 : 10 * i);
    end
end
cnt = 1;
for i = 1 : n
    for j = i : n
        Ym(cnt, :) = Ym2(j,:,i);
        cnt = cnt + 1;
    end
end

end

function A = A_matrix(V)
w = V(1:3);
v = V(4:6);
Aw = [zeros(3,1), -so_w(v), diag(w), [w(2), w(3), 0; w(1), 0, w(3); 0, w(1), w(2)]];
Av = [v, so_w(w), zeros(3,6)];
A = [Aw; Av];
end

