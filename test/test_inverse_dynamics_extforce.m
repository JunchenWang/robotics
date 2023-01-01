% robot = read_dynamics_file('F:\robotics\urdf\iiwa7\dynamics.txt');
% robot = convert_robot_tree(importrobot('E:\data\URDF\iiwa7\iiwa7.urdf'));
% robot = convert_robot_tree(importrobot('E:\data\URDF\ur_5e-calibrated\ur_description\urdf\ur5e-A302.urdf'));
robot = convert_robot_tree(importrobot('abbIrb120.urdf'));
[R, ~] = qr(rand(3));
robot.TCP = [R, [0, 0, 0.213]'; 0 0 0 1];
n = robot.dof;
q = 2 * pi * rand(1,n) - pi;
qd = 15 * (rand(1,n) - 0.5);
qdd = 10 * (rand(1,n) - 0.5);
extForce = 20 * (rand(6, n) - 0.5);
% extforce(:,:,1:6) = 0;
J = jacobian_matrix_all(robot, q);
ext_torque = zeros(n,1);
for i = 1 : n
    if i < n
        Tbc = [eye(3), robot.com(i,:)'; 0 0 0 1];
        extf = adjoint_T(Tbc)' * extForce(:,i);
    else
        extf = extForce(:,i);
    end
    ext_torque = ext_torque + J(:,:,i)'*extf;
end
[M, C, G] = mass_c_g_matrix(robot, q, qd);
tao1 = M * qdd' + C * qd' + G - ext_torque;
tao2 = inverse_dynamics_extforce(robot, q, qd, qdd, extForce);
disp(tao1);
disp(tao2);
disp(tao1 - tao2);