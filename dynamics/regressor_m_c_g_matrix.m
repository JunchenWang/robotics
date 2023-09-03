function [YTr, Mq, C, g, Jb, dJb, dMq, dTcp, Tcp] = regressor_m_c_g_matrix(robot, q, qd, a, v, r)
% paper: On the closed form computation of the dynamic matrices and their differentiations
% 注意，C矩阵和m_c_g_matrix不同，但都满足 dM = C + C'
n = robot.dof;
mass = robot.mass;
inertia = robot.inertia;
com = robot.com;
A = robot.A;
M = robot.M;
ME = robot.ME * robot.TCP;
J = zeros(6, n, n);
dJ = zeros(6, n, n);
gravity = [0;0;0;-robot.gravity'];
YTr = zeros(10 * n, 1);
% Y = zeros(n, 1);
C = zeros(n,n);
Mq = zeros(n, n);
g = zeros(n, 1);
dMq = zeros(n, n);
for i = 1 : n
    pik = [mass(i), mass(i) * com(i,:), inertia(1,1,i), inertia(2,2,i), inertia(3,3,i),...
                inertia(1,2,i), inertia(1,3,i), inertia(2,3,i)];
    Gk = spatial_inertia_matrix(inertia(:,:,i),mass(i), com(i,:));
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
    Mq = Mq + Jk' * Gk * Jk;
    dMq = dMq + dJk' * Gk * Jk + Jk' * Gk * dJk;
    C = C + Jk' * ((Gk * adk - adk' * Gk) * Jk + Gk * dJk);
    g = g + Jk' * Gk * Adk * gravity;
    alpha = Jk * a + adk * Jk * v + dJk * v +  rk;
    YTr((i - 1) * 10 + 1 : 10 * i) = (A_matrix(alpha) - adk' * A_matrix(Jk * v))' * Jk * r;
    % Y = Y + Jk' * (A_matrix(alpha) - adk' * A_matrix(Jk * v)) * pik';
    if i == n
        Tcp = tform_inv(T) * ME;
        dTcp = -tform_inv(T) * dT * tform_inv(T) * ME;
        adTb = adjoint_T(tform_inv(ME));
        Jb = adTb * Jk;
        dJb = adTb * dJk;
    end
end
% 
% [Mq2, C2, g2, Jb2, dJb2, dMq2, dTcp2, Tcp2] = m_c_g_matrix(robot, q, qd);
% 
% disp(norm(Y - Mq * a - C * v - g, 'fro'));
% disp(norm(Mq2 - Mq, 'fro'));
% disp(norm(C2 * qd - C * qd));
% disp(norm(g2 - g));
% disp(norm(Jb - Jb2, 'fro'));
% disp(norm(dJb - dJb2, 'fro'));
% disp(norm(dMq2 - dMq, 'fro'));
% disp(norm(dTcp2 - dTcp, 'fro'));
% disp(norm(Tcp2 - Tcp, 'fro'));
% disp(norm(dMq - C - C', 'fro'));
% disp(norm(dMq2 - C2 - C2', 'fro'));
% disp(norm(C - C2, 'fro'));
end

function A = A_matrix(V)
w = V(1:3);
v = V(4:6);
Aw = [zeros(3,1), -so_w(v), diag(w), [w(2), w(3), 0; w(1), 0, w(3); 0, w(1), w(2)]];
Av = [v, so_w(w), zeros(3,6)];
A = [Aw; Av];
end

