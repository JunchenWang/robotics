function gtao = gravity_torque(robot, q)
n = robot.dof;
gtao = inverse_dynamics(robot, q, zeros(n, 1), zeros(n, 1), zeros(6, 1));
