function test_task_space_dynamics
robot = convert_robot_tree(importrobot('URDF\iiwa7\iiwa7.urdf'));
q = -pi + 2 * pi * rand(1,7);
delta = 1e-8;
qd = -pi + 2 * pi * rand(1,7);
qdd = -pi + 2 * pi * rand(1,7);
Fext = rand(6, 7);
fext = -Fext(:,end);
ext_torque = get_ext_torque(robot, q, Fext);
[M, C, g, J, dJ, dM, dT, T] = m_c_g_matrix(robot, q, qd);
[M2, C2, g2, J2, dJ2, dM2, dT2, T2] = m_c_g_matrix(robot, q + delta * qd, qd + delta * qdd);
Z = null_z(J);
Z2 = null_z(J2);
dZ = derivative_null_z(J, dJ);
dZ2 = (Z2 - Z) / delta;
disp(norm(dZ - dZ2) / norm(dZ));

pinvZ = pinv_Z(Z, M);
pinvZ2 = pinv_Z(Z2, M2);
dpinvZ = derivative_pinv_Z(Z, M, dZ, dM);
dpinvZ2 = (pinvZ2 - pinvZ) / delta;
disp(norm(dpinvZ - dpinvZ2) / norm(dpinvZ));


pinvJ = pinv_J(J, M);
pinvJ2 = pinv_J(J2, M2);
dpinvJ = derivative_pinv_J(J, M, dJ, dM);
dpinvJ2 = (pinvJ2 - pinvJ) / delta;
disp(norm(dpinvJ - dpinvJ2) / norm(dpinvZ));

Ax = A_x(J, M);
Mux = Mu_x(J, M, dJ, C);
Av = A_v(Z, M);
Muv = Mu_v(Z, M, dZ, dM, C);
Muxv = Mu_xv(J, M, dJ, Z, C);
Muvx = Mu_vx(J, M, dZ, dM, Z, C);
dAx = -A_x_x(J, M, dJ * (M \ J') -J * (M \ dM) * (M \ J') + J * (M \ dJ')) * Ax;
dAv = dZ'*M*Z + Z' * dM * Z + Z' * M * dZ;

disp(norm(Muxv + Muvx'));
tem1 = dAx - 2 * Mux;
tem2 = dAv - 2 * Muv;
disp(norm(tem1 + tem1'));
disp(norm(tem2 + tem2'));

x = rand(6,1);
disp(norm(pinvJ * x - pinv_J_x(J,M,x)));
disp(norm(pinvJ' * q' - pinv_JT_x(J, M, q')));
