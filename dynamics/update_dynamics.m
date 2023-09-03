function robot = update_dynamics(robot, param)
n = robot.dof;
for i = 1 : n
    pi = param(10 * (i -1) + 1 : 10 * i);
    robot.mass(i) = pi(1);
    robot.com(i,:) = pi(2:4) / pi(1);
    robot.inertia(1,1,i) = pi(5);
    robot.inertia(2,2,i) = pi(6);
    robot.inertia(3,3,i) = pi(7);
    robot.inertia(1,2,i) = pi(8);
    robot.inertia(2,1,i) = pi(8);
    robot.inertia(1,3,i) = pi(9);
    robot.inertia(3,1,i) = pi(9);
    robot.inertia(2,3,i) = pi(10);
    robot.inertia(3,2,i) = pi(10);
end