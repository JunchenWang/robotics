function test_dynamics
robot = importrobot('E:\data\URDF\iiwa7\iiwa7.urdf');
% robot = importrobot('E:\data\URDF\ur_5e-calibrated\ur_description\urdf\ur5e-A302.urdf');
% robot = importrobot('kinovaGen3V12.urdf');
showdetails(robot);
flange = robot.Bodies{robot.NumBodies}.Name;
% flange = 'Bracelet_Link';
my_robot = convert_robot_tree(robot,flange);
[R, ~] = qr(rand(3));
my_robot.TCP = [R, [0, 0, 0.213]'; 0 0 0 1];
robot.Gravity = [0, 0, -9.8];
robot.DataFormat = 'row';
n = my_robot.dof;
q = rand(1,n) * 2 * pi - pi;
qd = rand(1,n) * 2 * pi - pi;
qdd = rand(1,n) * 2 * pi - pi;
gt1 = velocityProduct(robot, q, qd);
gt2 = velocity_torque(my_robot, q, qd);
disp(norm(gt1'-gt2));
F_ME = 100 * rand(6,1);
fext = externalForce(robot,flange,-F_ME,q);
F_ME = adjoint_T(my_robot.TCP)'*F_ME;
extf = zeros(6, n);
extf(:,end) = -F_ME;
tao = inverse_dynamics(my_robot, q, qd, qdd, F_ME);
tao3 = inverse_dynamics_extforce(my_robot, q, qd, qdd, extf);
qdd_x1 = forward_dynamics(my_robot, q, qd, tao, F_ME);
tao2 = inverseDynamics(robot, q, qd, qdd, fext);
qdd_x2 = forwardDynamics(robot, q, qd, tao2, fext);
disp(norm(tao'-tao2));
disp(norm(tao3'-tao2));
disp(norm(qdd_x1'-qdd_x2));

[Jb, T] = jacobian_matrix(my_robot, q);
Jb2 = geometricJacobian(robot, q, flange);
T2 = getTransform(robot, q, flange);
Jb2 = adjoint_T(tform_inv([T2(1:3,1:3), zeros(3,1); 0 0 0 1])) * Jb2;
T2 = T2 * my_robot.TCP;
Jb2 = adjoint_T(tform_inv(my_robot.TCP)) * Jb2;
disp(norm(T2-T));
disp(norm(Jb-Jb2));
% disp(T0);


