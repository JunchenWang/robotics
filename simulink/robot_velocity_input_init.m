rigid_body = importrobot('urdf\ur_description\urdf\ur5e-bph.urdf');
robot = get_robot_chain(rigid_body, 'tool0','world');
init_pos = deg2rad([0, -130, -100, -40, 90, 0]);