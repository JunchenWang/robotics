function M = CMatrix(robot, q, qd)
m = robot.m;
L = robot.L;
M = 0.5*m*L^2 * [-sin(q(2))*qd(2), -sin(q(2))*(qd(1)+qd(2)); sin(q(2))*qd(1), 0];
