function tao = inverse_dynamics_extforce(robot, q, qd, qdd, extForce)
% extForce is nx6 matrix, the wrench imposed by the envrionment to the
% robot's bodies
% extForce(:,end) is applied to TCP and others are applied to link frames
mass = robot.mass;
inertia = robot.inertia;
A = robot.A;
M = robot.M;
ME = robot.ME * robot.TCP;
com = robot.com;
n = robot.dof;
nu0 = zeros(6, 1);
dnu0 = [0, 0, 0, -robot.gravity]';
nu = zeros(6, n);
dnu = zeros(6, n);
tao = zeros(n, 1);
for i = 1 : n
    T = exp_twist(-A(i,:)*q(i))*tform_inv(M(:,:,i));
    Map = adjoint_T(T);
    nu(:, i) = Map * nu0 + A(i,:)'*qd(i);
    nu0 = nu(:, i);
    dnu(:, i) = Map * dnu0 + adjoint_twist(nu0') * A(i,:)' * qd(i) + A(i,:)'*qdd(i);
    dnu0 = dnu(:, i);
end

T = eye(4); 
F = zeros(6,1);
for i = n : -1 : 1
     if i < n
        Tbc = [eye(3), com(i,:)'; 0 0 0 1];
        extf = adjoint_T(Tbc)' * extForce(:,i);
     else
        extf = adjoint_T(tform_inv(ME))' * extForce(:,i);
     end
    G = [inertia(:,:,i), zeros(3);zeros(3), mass(i) * eye(3)];
    F = adjoint_T(T)'* F + G * dnu(:,i) - adjoint_twist(nu(:,i)')'*(G*nu(:,i)) - extf;
    tao(i) = F'*A(i,:)';
    T =  exp_twist(-A(i,:)*q(i))*tform_inv(M(:,:,i));
end