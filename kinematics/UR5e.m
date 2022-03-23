function robot = UR5e
% no dynamic data
a = [0.0, -0.425, -0.3922, 0.0, 0.0, 0.0];
d = [0.1625, 0.0, 0.0, 0.1333, 0.0997, 0.0996];
alpha = [1.570796327, 0, 0, 1.570796327, -1.570796327, 0];
delta_theta = [ -1.30565669838322851e-08, -0.588159115180872494, 0.826314252559053, -0.238165200723693349, -1.07362188606074938e-06, -2.93069160793579808e-07];
delta_a = [ 0.00018476643006531352, 0.0711540949987528037, 0.010821982777405037, -7.78445612723795073e-05, 3.66390856709752097e-05, 0];
delta_d = [ 0.000255399609815354367, -207.362134131258699, 220.009566418569648, -12.6472407036267196, 4.44834855784365102e-05, 5.03569516369511971e-05];
delta_alpha = [ -0.000158190453372508699, 0.00113799018245355158, 0.00732065515331056967, -0.000895020392787104768, 0.000581185017213048383, 0];
dh_table = [alpha(6) + delta_alpha(6), a(6) + delta_a(6), d(1) + delta_d(1), delta_theta(1);
           alpha(1) + delta_alpha(1), a(1) + delta_a(1), d(2) + delta_d(2), delta_theta(2);
           alpha(2) + delta_alpha(2), a(2) + delta_a(2), d(3) + delta_d(3), delta_theta(3);
           alpha(3) + delta_alpha(3), a(3) + delta_a(3), d(4) + delta_d(4), delta_theta(4);
           alpha(4) + delta_alpha(4), a(4) + delta_a(4), d(5) + delta_d(5), delta_theta(5);
           alpha(5) + delta_alpha(5), a(5) + delta_a(5), d(6) + delta_d(6), delta_theta(6)];
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
robot.ME = eye(4);
robot.mass = [3.761, 8.058, 2.846, 1.37, 1.3, 0.365]';
robot.inertia = zeros(3,3,6);
robot.jtMechanics = zeros(6,3);
robot.gravity=[0 0 -9.8000];