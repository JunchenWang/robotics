function [Mq, C, g] = mass_c_g_matrix(robot, q, qd)
mass = robot.mass;
inertia = robot.inertia;
A = robot.A;
M = robot.M;
n = robot.dof;
J = zeros(6, n, n);
Mq = zeros(n, n);
dJ = zeros(6, n, n, n);
dMq = zeros(n, n, n);
C = zeros(n, n);
G = zeros(6, 6, n);
P = zeros(n, n);
g = zeros(n,1);
for i = n : -1 : 1
    G(:,:,i) = [inertia(:,:,i), zeros(3);zeros(3), mass(i) * eye(3)];
    T = eye(4);
    dT = zeros([4,4,n]);
    for j = i : -1: 1
        tform = exp_twist(-A(j,:) * q(j)) * tform_inv(M(:,:,j));
        Tdtform = T * se_twist(-A(j,:)) * tform;
        for k = 1 : i
            [dAdT, AdT] = derivative_adjoint_T(T, dT(:,:,k));
            dJ(:,j,i,k)= dAdT * A(j,:)';
            dT(:,:,k) = dT(:,:,k) * tform;
            if k == j
                dT(:,:,k) = dT(:,:,k) + Tdtform;
            end
        end
        J(:,j,i) = AdT * A(j,:)';
        T = T * tform;
    end
    for k = 1 : i
        drc = -dT(1:3,1:3,k)'*T(1:3,4) - T(1:3,1:3)'*dT(1:3,4,k);
        P(i, k) = -mass(i) * dot(robot.gravity, drc);
    end
    Mq = Mq + J(:,:,i)'*G(:,:,i)*J(:,:,i);
end
for i = 1 : n
    for j = 1 : n
        dMq(:,:,i) = dMq(:,:,i) + dJ(:,:,j,i)'*G(:,:,j)*J(:,:,j) + J(:,:,j)'*G(:,:,j)*dJ(:,:,j,i);
        g(i) = g(i) + P(j, i);
    end
end
for k = 1 : n
    for j = 1 : n
        for i = 1 : n
            C(k, j) = C(k, j) + 0.5 * (dMq(k,j,i) + dMq(k,i,j) - dMq(i,j,k))*qd(i);
        end
    end
end