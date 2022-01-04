function [Mq, J] = mass_matrix2(robot, q)
mass = robot.mass;
inertia = robot.inertia;
A = robot.A;
M = robot.M;
n = robot.dof;
J = zeros(6, n, n);
Mq = zeros(n, n);
for i = 1 : n
   G = [inertia(:,:,i), zeros(3);zeros(3), mass(i) * eye(3)];
   T =  exp_twist(-A(i,:) * q(i)) * tform_inv(M(:,:,i));
   J(:, i, i)  = A(i,:)';
   AdT = adjoint_T(T);
   for j = 1 : i - 1
       A(j,:) = AdT * A(j,:)';
       J(:, j, i)  = A(j,:)';
   end
   Mq = Mq + J(:,:,i)'*G*J(:,:,i);
end