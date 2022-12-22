function tao = inverse_dynamics_extforce(robot, q, qd, qdd, extForce)
% extForce is nx6 matrix, the wrench imposed by the envrionment to the robot's end-effector
mass = robot.mass;
inertia = robot.inertia;
A = robot.A;
M = robot.M;
ME = robot.ME;
%jtMechanics = robot.jtMechanics;
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
extForce(:,end) = adjoint_T(tform_inv(ME))' * extForce(:,end);
T = eye(4); 
F = zeros(6,1);
for i = n : -1 : 1
    G = [inertia(:,:,i), zeros(3);zeros(3), mass(i) * eye(3)];
    F = adjoint_T(T)'* F + G * dnu(:,i) - adjoint_twist(nu(:,i)')'*(G*nu(:,i)) - extForce(:, i);
    tao(i) = F'*A(i,:)'; %+ jtMechanics(i, 1) * qd(i) + jtMechanics(i, 2) * (q(i) - jtMechanics(i, 3)); % consider joint mechanics
    T =  exp_twist(-A(i,:)*q(i))*tform_inv(M(:,:,i));
end