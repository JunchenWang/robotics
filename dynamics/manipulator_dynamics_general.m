function yd = manipulator_dynamics_general(robot, controller, fext, t, y)
% fext is a force function applied to the robot
% fext(t,y)(:,end) is applied to TCP and others are applied to link frames
n = robot.dof;
q = y(1:n);
qd = y(n + 1 : 2 * n);
hqqd = gravity_velocity_torque(robot, q, qd);
Fext = fext(t, y);
ext_torque = get_ext_torque(robot, q, Fext);
[tau, td] = controller(t, y, Fext);
sz = length(td);
yd = zeros(2 * n + sz, 1);
yd(1:n) = qd;
yd(n + 1 : 2 * n) = tau - hqqd + ext_torque;
yd(2 * n + 1 : end) = td;
end