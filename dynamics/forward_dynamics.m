function qdd = forward_dynamics(robot, q, qd, tao, F_ME)
ME = robot.ME * robot.TCP;
n = robot.dof;
[Mq, J] = mass_matrix(robot, q);
Jb = adjoint_T(tform_inv(ME)) * J(:,:,n);
hqqd = inverse_dynamics(robot, q, qd, zeros(n, 1), zeros(6, 1));
b = tao - hqqd - Jb' * F_ME;
qdd = Mq \ b;
% M*qdd = b
% disp(Jb);
