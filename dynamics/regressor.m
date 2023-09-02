function Y = regressor(robot, q, qd, a, v)
n = robot.dof;
mass = robot.mass;
inertia = robot.inertia;
com = robot.com;
A = robot.A;
M = robot.M;
J = zeros(6, n, n);
dJ = zeros(6, n, n);
g = -robot.gravity';
Y = zeros(n,1);
for i = 1 : n
    ppi = [mass(i), mass(i) * com(i,:), inertia(1,1,i), inertia(2,2,i), inertia(3,3,i),...
                inertia(1,2,i), inertia(1,3,i), inertia(2,3,i)];
    T = eye(4);
    dT = zeros(4,4);% T对时间的导数
    for j = i : -1: 1
        [dAdT, AdT] = derivative_adjoint_T(T, dT);
        J(:,j,i) = AdT * A(j,:)';
        dJ(:,j,i)= dAdT * A(j,:)';
        tform = exp_twist(-A(j,:) * q(j)) * tform_inv(M(:,:,j));
        dT = (dT  + T * se_twist(-A(j,:)) * qd(j)) * tform;
        T = T * tform;
    end
    gamma = adjoint_T(T) * [zeros(3,1);g];
    V = J(:,:,i) * qd;
    alpha = J(:,:,i) * a + adjoint_twist(V) * J(:,:,i) * v...
            + dJ(:,:,i) * v +  gamma;
    Y = Y + J(:,:,i)' * (A_matrix(alpha) - adjoint_twist(V)' * A_matrix(J(:,:,i) * v)) * ppi';
    % disp(G(:,:,i) * V - A_matrix(V) * ppi(:,i));
end

[Mq, C, G] = m_c_g_matrix(robot, q, qd);
disp(Y - Mq * a - C * v - G);
end

function A = A_matrix(V)
w = V(1:3);
v = V(4:6);
Aw = [zeros(3,1), -so_w(v), diag(w), [w(2), w(3), 0; w(1), 0, w(3); 0, w(1), w(2)]];
Av = [v, so_w(w), zeros(3,6)];
A = [Aw; Av];
end

