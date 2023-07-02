function yd = manipulator_dynamics_observer(robot, controller, fext, t, y)
% fext is a force function applied to the robot
% fext(t,y)(:,end) is applied to TCP and others are applied to link frames
n = robot.dof;
yd = zeros(3*n,1);
q = y(1:n);
qd = y(n + 1 : 2 * n);
hqqd = gravity_velocity_torque(robot, q, qd);
yd(1:n) = qd;
ext_torque = get_ext_torque(robot, q, fext(t, y));
[tau, td] = controller(t, y);
yd(n + 1 : 2 * n) = tau - hqqd + ext_torque;
yd(2 * n + 1 : end) = td;