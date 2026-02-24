function robot = read_urdf(filename)

robot = get_robot_chain(importrobot(filename));

end