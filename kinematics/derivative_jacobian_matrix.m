function [dJb, Jb, dT, T] = derivative_jacobian_matrix(robot, q, qd)
A = robot.A;
M = robot.M;
ME = robot.ME;
n = robot.dof;
T = ME;
Jb = zeros(6, n);
dJb = zeros(6, n);
dT = zeros(4, 4);
for i = n : -1 : 1
    [dInvT, invT] = derivative_tform_inv(T, dT);
    [dAdT, adT] = derivative_adjoint_T(invT, dInvT);
    Jb(:,i) = adT * A(i,:)';
    dJb(:, i) = dAdT * A(i,:)';
    tform = M(:,:,i) * exp_twist(A(i,:) * q(i));
    dT = tform * (se_twist(A(i,:)) * T * qd(i) + dT);
    T = tform * T;
end