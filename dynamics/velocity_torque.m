function vtao = velocity_torque(robot, q, qd)
gvtao = gravity_velocity_torque(robot, q, qd);
vtao = gvtao - gravity_torque(robot, q);
