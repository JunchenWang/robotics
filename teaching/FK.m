function [x, y, theta] = FK(robot, q)
L = robot.L;
x = -L * (sin(q(1)) + sin(q(1) + q(2)));
y = L * (cos(q(1)) + cos(q(1) + q(2)));
theta = q(1) + q(2);