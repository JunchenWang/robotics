function G = GMatrix(robot, q)
m = robot.m;
L = robot.L;
g = robot.g;
G = -0.5*m*g*L*[3*sin(q(1))+sin(q(1) + q(2));sin(q(1) + q(2))];
