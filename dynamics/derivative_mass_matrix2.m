function [dMq, Mq, dJ, J] = derivative_mass_matrix2(robot, q, qd)
mass = robot.mass;
inertia = robot.inertia;
A = robot.A;
M = robot.M;
n = robot.dof;
J = zeros(6, n, n);
dJ = zeros(6, n, n);
Mq = zeros(n, n);
dMq = zeros(n, n);
dA = zeros(size(A));
for i = 1 : n
   G = [inertia(:,:,i), zeros(3);zeros(3), mass(i) * eye(3)];
   T =  exp_twist(-A(i,:) * q(i)) * tform_inv(M(:,:,i));
   dT = se_twist(-A(i,:)) * qd(i) * T;
   J(:, i, i)  = A(i,:)';
   [dAdT, AdT] = derivative_adjoint_T(T, dT);
   for j = 1 : i - 1
       dA(j,:) = dAdT * A(j,:)' + AdT * dA(j, :)';
       A(j,:) = AdT * A(j,:)';
       J(:, j, i)  = A(j,:)';
       dJ(:, j, i)  = dA(j,:)';
   end
   Mq = Mq + J(:,:,i)'*G*J(:,:,i);
   dMq = dMq + dJ(:,:,i)'*G*J(:,:,i) + J(:,:,i)'*G*dJ(:,:,i);
end