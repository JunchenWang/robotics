function [dMq, Mq, dJ, J] = derivative_mass_matrix(robot, q, qd)
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
% 这个方法计算dJ不对，目前还没找到原因，derivative_mass_matrix2 的结果是对的
for i = 1 : n
   G = [inertia(:,:,i), zeros(3);zeros(3), mass(i) * eye(3)];
   T =  exp_twist(-A(i,:) * q(i)) * tform_inv(M(:,:,i));
   dT = se_twist(-A(i,:)) * qd(i) * T;
   J(:, i, i)  = A(i,:)';
   [dAdT, AdT] = derivative_adjoint_T(T, dT);
   for j = 1 : i - 1
       A(j,:) = AdT * A(j,:)';
       dA(j,:) = dAdT * A(j,:)' + AdT * dA(j, :)';
       J(:, j, i)  = A(j,:)';
       dJ(:, j, i)  = dA(j,:)';
   end
   Mq = Mq + J(:,:,i)'*G*J(:,:,i);
   dMq = dMq + dJ(:,:,i)'*G*J(:,:,i) + J(:,:,i)'*G*dJ(:,:,i);
end