function C = c_matrix(Mq, q, qd)
% Jb, dJb expressed in TCP
mass = robot.mass;
com = robot.com;
inertia = robot.inertia;
A = robot.A;
M = robot.M;
ME = robot.ME * robot.TCP;
n = robot.dof;
J = zeros(6, n, n);
dJ = zeros(6, n, n);
Mq = zeros(n, n);
dMq = zeros(n, n);
pdJ = zeros(6, n, n, n);%关于关节角的偏导数
pdMq = zeros(n, n, n);%关于关节角的偏导数
C = zeros(n, n);
G = zeros(6, 6, n);
P = zeros(n, n);
g = zeros(n,1);
Jb = zeros(6, n);
dJb = zeros(6, n);
Tcp = eye(4);
dTcp = zeros(4);
for i = n : -1 : 1
    G(:,:,i) = spatial_inertia_matrix(inertia(:,:,i),mass(i), com(i,:));
    T = eye(4);% 连杆i的坐标系的逆
    dT = zeros(4,4);% T对时间的导数
    pdT = zeros(4,4,n);% T对关节变量的偏导数
    for j = i : -1: 1
        [dAdT, AdT] = derivative_adjoint_T(T, dT);
        J(:,j,i) = AdT * A(j,:)';
        dJ(:,j,i)= dAdT * A(j,:)';
        tform = exp_twist(-A(j,:) * q(j)) * tform_inv(M(:,:,j));
        % tform 对qj的偏导数
        pdtform = se_twist(-A(j,:)) * tform;
        % 雅克比矩阵对关节变量的偏导数
        pdT(:,:,j) = T * pdtform;
        for k = j + 1 : i
            pdAdT = derivative_adjoint_T(T, pdT(:,:,k));
            pdJ(:,j,i,k)= pdAdT * A(j,:)';
            pdT(:,:,k) = pdT(:,:,k) * tform;
            % if k == j
            %     pdT(:,:,k) = pdT(:,:,k) + T * pdtform;
            % end
        end
        dT = (dT  + T * se_twist(-A(j,:)) * qd(j)) * tform;
        T = T * tform;
    end
    Mq = Mq + J(:,:,i)'*G(:,:,i)*J(:,:,i);
    dMq = dMq + dJ(:,:,i)'*G(:,:,i)*J(:,:,i) + J(:,:,i)'*G(:,:,i)*dJ(:,:,i);
    for k = 1 : i
        drc = -pdT(1:3,1:3,k)'*(T(1:3,4) - com(i,:)') - T(1:3,1:3)'*pdT(1:3,4,k); %T的逆是第i个连杆的坐标系矩阵
        P(i, k) = -mass(i) * dot(robot.gravity, drc); %注意负号
    end
    if i == n
        Tcp = tform_inv(T) * ME;
        dTcp = -tform_inv(T) * dT * tform_inv(T) * ME;
        adTb = adjoint_T(tform_inv(ME));
        Jb = adTb * J(:,:,n);
        dJb = adTb * dJ(:,:,n);
    end
end
% g matrix
for i = 1 : n
    for j = 1 : n
        pdMq(:,:,i) = pdMq(:,:,i) + pdJ(:,:,j,i)'*G(:,:,j)*J(:,:,j) + J(:,:,j)'*G(:,:,j)*pdJ(:,:,j,i);
        g(i) = g(i) + P(j, i);
    end
end
% C matrix
for k = 1 : n
    for j = 1 : n
        for i = 1 : n
            C(k, j) = C(k, j) + 0.5 * (pdMq(k,j,i) + pdMq(k,i,j) - pdMq(i,j,k))*qd(i);
        end
    end
end