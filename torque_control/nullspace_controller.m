function tao = nullspace_controller(robot, t, y, Td, kesai, Dn, Kn)
% error in desired frame
n = robot.dof;
q = y(1:n);
qd = y(n + 1 : 2 * n);
q0 = inverse_kin_kuka_robot_kesai_near(robot, Td, kesai, q);
tao = -Dn * qd - Kn * (q - q0');
[Jb, T] = jacobian_matrix(robot, q);
M = mass_matrix(robot, q);
R = T(1:3,1:3);
Rd = Td(1:3, 1:3);
re = logR(Rd'*R)';
Ar = w_dr_A(re);
Jw = Jb(1:3,:);
Jv = Jb(4:6,:);
Jx = [Ar \ Jw; Rd'*R*Jv];
tao = tao - Jx' * ((Jx * (M \ Jx'))  \ (Jx * (M \ tao)));