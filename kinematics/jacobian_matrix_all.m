function J = jacobian_matrix_all(robot, q)
A = robot.A;
M = robot.M;
n = robot.dof;
com = robot.com;
J = zeros(6, n, n);
ME = robot.ME * robot.TCP;
Tb = eye(4);
for i = n : -1 : 1
   T = eye(4);
   for j = i : -1: 1
       J(:,j,i) = adjoint_T(T) * A(j,:)';
       T = T * exp_twist(-A(j,:) * q(j)) * tform_inv(M(:,:,j));
   end
   if i < n
       Tb(1:3,4) = com(i,:);
       J(:,:,i) = adjoint_T(Tb) * J(:,:,i);
   end
end
J(:,:,end) = adjoint_T(tform_inv(ME)) * J(:,:,end);