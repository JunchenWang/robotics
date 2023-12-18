function M = mass_matrix(robot, q)
m = robot.m;
L = robot.L;
M = m*L^2 * [5/3 + cos(q(2)), 1/3+0.5*cos(q(2)); 1/3+0.5*cos(q(2)), 1/3];
end
