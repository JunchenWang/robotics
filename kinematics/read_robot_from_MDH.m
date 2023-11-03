function robot = read_robot_from_MDH(dh_table)
% dh_table: alpha_i-1, a_i-1, d_i, theta_i

n = size(dh_table, 1);
robot.dof = n;
M = zeros(4,4,n);
A = zeros(n, 6);
for i = 1 : 6
    M(:,:,i) = rttr(dh_table(i,:));
    A(i,:) = [0, 0, 1, 0, 0, 0];
end
robot.M = M;
robot.A = A;
robot.ME = eye(4);
robot.mass = zeros(1,n);
robot.inertia = zeros(3,3,6);
robot.com = zeros(6,3);
robot.gravity=[0 0 -9.8000];
robot.TCP = eye(4);