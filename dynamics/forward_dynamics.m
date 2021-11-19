function qdd = forward_dynamics(robot, q, qd, tao, F_ME)
mass = robot.mass;
inertia = robot.inertia;
A = robot.A;
M = robot.M;
ME = robot.ME;
n = robot.dof;
J = zeros(6, n, n);
Mq = zeros(n, n);
for i = 1 : n
   G = [inertia(:,:,i), zeros(3);zeros(3), mass(i) * eye(3)];
   T = M(:,:,i) * exp_twist(A(i,:) * q(i));
   J(:, i, i)  = A(i,:)';
   for j = 1 : i - 1
       A(j,:) = adjoint_T(tform_inv(T)) * A(j,:)';
       J(:, j, i)  = A(j,:)';
   end
   Mq = Mq + J(:,:,i)'*G*J(:,:,i);
end
Jb = adjoint_T(tform_inv(ME)) * J(:,:,n);
hqqd = inverse_dynamics(robot, q, qd, zeros(n, 1), zeros(6, 1));
b = tao - hqqd - Jb' * F_ME;
qdd = Mq \ b;
% M*qdd = b
% disp(Jb);
