function dynamics
iiwa = read_dynamics_file('F:\MICR\MICSys\dynamics.txt');
% TCP = eye(4);
% TCP(1:3,4) = [0, 0, 0.1];
% centerofmass = [0, 0, 0]';
% load = 0;
% R = eye(3);
% I = zeros(3);
% [mass, inertia, A, M, ME] = attach_tool(mass, inertia, A, M, ME, TCP, centerofmass, load, R, I);
robot = importrobot('E:\data\URDF\iiwa7\iiwa7.urdf');
robot.Gravity = [0, 0, -9.8];
robot.DataFormat = 'row';
% q = [10, 20, 30, 40, 50, 60, 70] / 180 * pi;
% q = zeros(7,1);
% qd = zeros(7,1);
% qdd = zeros(7,1);
% F_ME = zeros(6,1);

q = rand(1,7) * 2 * pi - pi;
qd = rand(1,7) * 2 * pi - pi;
qdd = rand(1,7) * 2 * pi - pi;
gt1 = velocityProduct(robot, q, qd);
gt2 = velocity_torque(iiwa, q, qd);
disp(norm(gt1'-gt2));
F_ME = 100 * rand(6,1);
fext = externalForce(robot,'iiwa_link_ee_kuka',-F_ME,q);
tao = inverse_dynamics(iiwa, q, qd, qdd, F_ME);
qdd_x1 = forward_dynamics(iiwa, q, qd, tao, F_ME);
tao2 = inverseDynamics(robot, q, qd, qdd, fext);
qdd_x2 = forwardDynamics(robot, q, qd, tao2, fext);
disp(norm(tao'-tao2));
disp(norm(qdd_x1'-qdd_x2));

[Jb, T] = jacobian_matrix(iiwa, q);
[B, T0] = jacobian_matrix(iiwa, zeros(1,7));
Jb2 = jacobian(B', q, 'b');
T2 = forward_kin_kuka(q);
T2(1:3,4) = T2(1:3,4) / 1000;
disp(norm(T2-T));
disp(norm(Jb-Jb2));


