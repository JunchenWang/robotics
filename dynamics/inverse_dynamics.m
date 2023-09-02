function tao = inverse_dynamics(robot, q, qd, qdd, F_ME)
% F_ME is the wrench imposed to the envrionment by the robot's TCP
% note that the difference with inverse_dynamics_extforce
mass = robot.mass;
inertia = robot.inertia;
com = robot.com;
A = robot.A;
M = robot.M;
ME = robot.ME * robot.TCP;
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
F = F_ME;
T = tform_inv(ME);
for i = n : -1 : 1
    G = spatial_inertia_matrix(inertia(:,:,i),mass(i), com(i,:));
    F = adjoint_T(T)'* F + G * dnu(:,i) - adjoint_twist(nu(:,i)')'*(G*nu(:,i));
    tao(i) = F'*A(i,:)'; %+ jtMechanics(i, 1) * qd(i) + jtMechanics(i, 2) * (q(i) - jtMechanics(i, 3)); % consider joint mechanics
    T =  exp_twist(-A(i,:)*q(i))*tform_inv(M(:,:,i));
end