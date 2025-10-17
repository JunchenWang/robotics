function dM = d_mass_matrix(robot, q, qd)
m = robot.m;
L = robot.L;
dM = m*L^2 * [-sin(q(2))*qd(2),   -0.5*sin(q(2))*qd(2);-0.5*sin(q(2))*qd(2), 0];
end
