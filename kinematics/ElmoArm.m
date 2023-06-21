function robot = ElmoArm
% no dynamic data
alpha = [pi, pi/2, pi, pi, pi/2, pi/2];
a = [0, 0, 0.385, 0.315, 0, 0];
d = [-0.136,  -0.1245, 0, 0, 0.095, -0.0895];
delta_theta = [pi, -pi/2, pi/2, pi/2, pi/2, pi];
delta_a = [ 0, 0, 0, 0, 0, 0];
delta_d = [0, 0, 0, 0, 0, 0];
delta_alpha = [0, 0, 0, 0, 0, 0];
dh_table = [alpha(1) + delta_alpha(1), a(1) + delta_a(1), d(1) + delta_d(1), delta_theta(1);
           alpha(2) + delta_alpha(2), a(2) + delta_a(2), d(2) + delta_d(2), delta_theta(2);
           alpha(3) + delta_alpha(3), a(3) + delta_a(3), d(3) + delta_d(3), delta_theta(3);
           alpha(4) + delta_alpha(4), a(4) + delta_a(4), d(4) + delta_d(4), delta_theta(4);
           alpha(5) + delta_alpha(5), a(5) + delta_a(5), d(5) + delta_d(5), delta_theta(5);
           alpha(6) + delta_alpha(6), a(6) + delta_a(6), d(6) + delta_d(6), delta_theta(6)];
n = 6;
robot.dof = n;
M = zeros(4,4,n);
A = zeros(n, 6);
for i = 1 : 6
    M(:,:,i) = rttr(dh_table(i,:));
    A(i,:) = [0, 0, 1, 0, 0, 0];
end
robot.M = M;
robot.A = A;
robot.ME = [1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1];
robot.mass = [3.761, 8.058, 2.846, 1.37, 1.3, 0.365]';
robot.inertia = zeros(3,3,6);
robot.com = zeros(6,3);
robot.gravity=[0 0 -9.8000];
robot.TCP = eye(4);