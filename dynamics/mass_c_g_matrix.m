function [Mq, C, g, Jb, dJb, Tcp] = mass_c_g_matrix(robot, q, qd)
% Jb, dJb expressed in TCP
mass = robot.mass;
inertia = robot.inertia;
A = robot.A;
M = robot.M;
ME = robot.ME * robot.TCP;
com = robot.com;
n = robot.dof;
J = zeros(6, n, n);
Mq = zeros(n, n);
pdJ = zeros(6, n, n, n);%关于关节角的偏导数
pdMq = zeros(n, n, n);%关于关节角的偏导数
C = zeros(n, n);
G = zeros(6, 6, n);
P = zeros(n, n);
g = zeros(n,1);
Jb = zeros(6, n);
dJb = zeros(6, n);
for i = n : -1 : 1
    G(:,:,i) = spatial_inertia_matrix(inertia(:,:,i),mass(i), com(i,:));
    T = eye(4);
    dT = zeros([4,4,n]);
    for j = i : -1: 1
        tform = exp_twist(-A(j,:) * q(j)) * tform_inv(M(:,:,j));
        Tdtform = T * se_twist(-A(j,:)) * tform;
        for k = 1 : i
            [dAdT, AdT] = derivative_adjoint_T(T, dT(:,:,k));
            pdJ(:,j,i,k)= dAdT * A(j,:)';
            dT(:,:,k) = dT(:,:,k) * tform;
            if k == j
                dT(:,:,k) = dT(:,:,k) + Tdtform;
            end
        end
        J(:,j,i) = AdT * A(j,:)';
        T = T * tform;
    end
    for k = 1 : i
        drc = -dT(1:3,1:3,k)'*(T(1:3,4) - com(i,:)') - T(1:3,1:3)'*dT(1:3,4,k);
        P(i, k) = -mass(i) * dot(robot.gravity, drc);
    end
    Mq = Mq + J(:,:,i)'*G(:,:,i)*J(:,:,i);
    if i == n
        Tcp = tform_inv(T) * ME;
        Jb = J(:,:,n);
    end
end

for i = 1 : n
    dJb = dJb + pdJ(:,:,n,i) * qd(i);
    for j = 1 : n
        pdMq(:,:,i) = pdMq(:,:,i) + pdJ(:,:,j,i)'*G(:,:,j)*J(:,:,j) + J(:,:,j)'*G(:,:,j)*pdJ(:,:,j,i);
        g(i) = g(i) + P(j, i);
    end
end
% jacobian and its derivatives
Tb = tform_inv(ME);
adTb = adjoint_T(Tb);
Jb = adTb * Jb;
dJb = adTb * dJb;
% C matrix
for k = 1 : n
    for j = 1 : n
        for i = 1 : n
            C(k, j) = C(k, j) + 0.5 * (pdMq(k,j,i) + pdMq(k,i,j) - pdMq(i,j,k))*qd(i);
        end
    end
end