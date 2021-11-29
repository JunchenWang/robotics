function tao = inverse_dynamics_task(robot, q, V, dV, F_ME)
% n = robot.dof;
Mq = mass_matrix(robot, q);
Jb = jacobian_matrix(robot, q);
qd = lsqminnorm(Jb, V);
[~, ~, dJb, ~] = derivative_jacobian_matrix(robot, q, qd);
Jb_inv = pinv(Jb);
%%
% Aq = Jb_inv' * Mq * Jb_inv;
% hqqd = gravity_velocity_torque(robot, q, qd);
% etaqv = Jb_inv' * hqqd - Aq*dJb*qd;
% F = Aq * dV + etaqv + F_ME;
%% 
% Aq =  Mq * Jb_inv;
% hqqd = gravity_velocity_torque(robot, q, qd);
% etaqv = hqqd - Aq*dJb*qd;
% tao = Aq * dV + etaqv + Jb' * F_ME;
%%
hqqd = gravity_velocity_torque(robot, q, qd);
tao = Mq * (Jb_inv * dV - Jb_inv * dJb * qd) + hqqd + Jb'*F_ME;