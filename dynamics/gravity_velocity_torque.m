function gvtao = gravity_velocity_torque(robot, q, qd)
n = robot.dof;
gvtao = inverse_dynamics(robot, q, qd, zeros(n, 1), zeros(6, 1));

