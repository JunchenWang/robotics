rigid_body = importrobot('urdf\ur_description\urdf\ur5e.urdf');
robot = get_robot_chain(rigid_body, 'tool0','world');
init_pos = deg2rad([0, -130, -100, -40, 90, 0]);
T_init = forward_kin_general(robot, init_pos);