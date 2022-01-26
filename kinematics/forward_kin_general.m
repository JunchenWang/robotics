function T = forward_kin_general(robot, q)
A = robot.A;
M = robot.M;
ME = robot.ME;
n = robot.dof;
T = ME;
for i = n : -1 : 1
    T = M(:,:,i) * exp_twist(A(i,:) * q(i)) * T;
end