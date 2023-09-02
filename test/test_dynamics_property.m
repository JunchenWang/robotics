function test_dynamics_property
robot = convert_robot_tree2(importrobot('URDF\iiwa7\iiwa7.urdf'));
q = -pi + 2 * pi * rand(1,7);
qd = -pi + 2 * pi * rand(1,7);
qdd = -pi + 2 * pi * rand(1,7);
Fext = rand(6, 7);
fext = -Fext(:,end);
ext_torque = get_ext_torque(robot, q, Fext);
% tic;
% [M2, C2, g2, J2, dJ2, dM2, dT2, T2] = m_c_g_matrix2(robot, q, qd);
% toc;
% tic;
[M, C, g, J, dJ, dM, dT, T] = m_c_g_matrix(robot, q, qd);
% toc;
% err = norm(M2 - M, 'fro') + norm(C2 - C, 'fro') +norm(g2 - g, 'fro') ...
%       + norm(dJ2 - dJ, 'fro')  + norm(dM2 - dM, 'fro')  + norm(dT2 - dT, 'fro')...
%       + norm(T2 - T, 'fro') ;
% disp(err);
[M3, C3, g3, J3, dJ3, T3] = mass_c_g_matrix(robot, q, qd);
[dM1, M1, dJ1, J1] = derivative_mass_matrix(robot,q,qd);
TCP = robot.ME * robot.TCP;
dJ1 = adjoint_T(tform_inv(TCP)) * dJ1(:,:,end);
J1 = adjoint_T(tform_inv(TCP)) * J1(:,:,end);
[dJ2, J2, dT2, T2] = derivative_jacobian_matrix(robot, q, qd);
% tao1 = inverse_dynamics(robot, q, qd, qdd, fext);
% tic;
tao2 = inverse_dynamics_fext(robot, q, qd, qdd, Fext);
% toc;
tao3 = M * qdd' + C * qd' + g - ext_torque;

disp(norm(tao2 - tao3));
disp(norm(M - M1, 'fro'));
disp(norm(dM - dM1, 'fro'));
disp(norm(J - J1, 'fro'));
disp(norm(dJ - dJ1, 'fro'));
disp(norm(J - J2, 'fro'));
disp(norm(dJ - dJ2, 'fro'));
disp(norm(T - T2, 'fro'));
disp(norm(dT - dT2, 'fro'));

disp(norm(M - M3, 'fro'));
disp(norm(C - C3, 'fro'));
disp(norm(g - g3, 'fro'));

disp(norm(J - J3, 'fro'));
disp(norm(dJ - dJ3, 'fro'));
disp(norm(T - T3, 'fro'));

disp(norm(dM - (C + C'), 'fro')); % dM = C + C'; dM - 2*C is skew-symmetric
% disp(norm(M));
