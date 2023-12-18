function q = two_link_IK(robot, x, y, cfg)
L = robot.L;
R2 = x^2+y^2;
if R2 > 4 * L^2
    error('no solution');
end
c2 = (R2/L^2 - 2) / 2;
q = zeros(1,2);
q(2) = cfg * acos(c2);
a = acos(R2 / (2 * sqrt(R2) * L));
q(1) = atan2(-x, y) - cfg * a;
end