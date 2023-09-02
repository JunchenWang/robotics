function [Mq, J] = mass_matrix(robot, q)
% J is expressed in link com, different with jacobian_matrix
mass = robot.mass;
inertia = robot.inertia;
A = robot.A;
M = robot.M;
n = robot.dof;
com = robot.com;
J = zeros(6, n, n);
Mq = zeros(n, n);
for i = n : -1 : 1
   G = spatial_inertia_matrix(inertia(:,:,i),mass(i), com(i,:));
   T = eye(4);
   for j = i : -1: 1
       J(:,j,i) = adjoint_T(T) * A(j,:)';
       T = T * exp_twist(-A(j,:) * q(j)) * tform_inv(M(:,:,j));
   end
   Mq = Mq + J(:,:,i)'*G*J(:,:,i);
end