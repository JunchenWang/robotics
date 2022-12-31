function [Jb, T] = jacobian_matrix(robot, q)
A = robot.A;
M = robot.M;
ME = robot.ME * robot.TCP;
n = robot.dof;
T = ME;
Jb = zeros(6, n);
for i = n : -1 : 1
    Jb(:,i) = adjoint_T(tform_inv(T)) * A(i,:)';
    T = M(:,:,i) * exp_twist(A(i,:) * q(i)) * T;
end