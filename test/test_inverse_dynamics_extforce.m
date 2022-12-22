robot = read_dynamics_file('F:\robotics\urdf\iiwa7\dynamics.txt');
q = 2 * pi * rand(1,7) - pi;
qd = 15 * (rand(1,7) - 0.5);
qdd = 10 * (rand(1,7) - 0.5);
extforce = 20 * (rand(6, 7) - 0.5);
% extforce(:,:,1:6) = 0;
J = jacobian_matrix_all(robot, q);
ext_torque = zeros(7,1);
for i = 1 : 7
    ext_torque = ext_torque + J(:,:,i)'*extforce(:,i);
end
[M, C, G] = mass_c_g_matrix(robot, q, qd);
tao1 = M * qdd' + C * qd' + G - ext_torque;
tao2 = inverse_dynamics_extforce(robot, q, qd, qdd, extforce);
disp(tao1);
disp(tao2);
tao1 - tao2