function [dMq, Mq, dJ, J] = derivative_mass_matrix(robot, q, qd)
% J is expressed in link com, different with jacobian_matrix
mass = robot.mass;
inertia = robot.inertia;
com = robot.com;
A = robot.A;
M = robot.M;
n = robot.dof;
J = zeros(6, n, n);
Mq = zeros(n, n);
dJ = zeros(6, n, n);
dMq = zeros(n, n);
for i = n : -1 : 1
   G = spatial_inertia_matrix(inertia(:,:,i),mass(i), com(i,:));
   T = eye(4);
   dT = zeros(4);
   for j = i : -1: 1
       [dAdT, AdT] = derivative_adjoint_T(T, dT);
       J(:,j,i) = AdT * A(j,:)';
       dJ(:,j,i)= dAdT * A(j,:)';
       tform = exp_twist(-A(j,:) * q(j)) * tform_inv(M(:,:,j));
       dT = (dT  + T * se_twist(-A(j,:)) * qd(j)) * tform;
       T = T * tform;
   end
   Mq = Mq + J(:,:,i)'*G*J(:,:,i);
   dMq = dMq + dJ(:,:,i)'*G*J(:,:,i) + J(:,:,i)'*G*dJ(:,:,i);
end