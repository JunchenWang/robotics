function yd = manipulator_dynamics_fext(robot, tao, fext, t, y)
% fext is a force function applied to the robot
% fext(t,y)(:,end) is applied to TCP and others are applied to link frames
n = robot.dof;
yd = zeros(2*n,1);
q = y(1:n);
qd = y(n + 1 : end);
hqqd = gravity_velocity_torque(robot, q, qd);
yd(1:n) = qd;
ext_torque = get_ext_torque(robot, q, fext(t, y));
yd(n + 1 : end) = tao(t, y) - hqqd + ext_torque;