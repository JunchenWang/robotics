function J = jacobian_matrix_all(robot, q)
A = robot.A;
M = robot.M;
n = robot.dof;
J = zeros(6, n, n);
for i = n : -1 : 1
   T = eye(4);
   for j = i : -1: 1
       J(:,j,i) = adjoint_T(T) * A(j,:)';
       T = T * exp_twist(-A(j,:) * q(j)) * tform_inv(M(:,:,j));
   end
end
J(:,:,end) = adjoint_T(tform_inv(robot.ME)) * J(:,:,end);