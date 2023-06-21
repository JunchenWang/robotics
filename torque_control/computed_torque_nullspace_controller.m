function tao = computed_torque_nullspace_controller(robot, t, y, Td, kesai, Dn, Kn, freq)
n = robot.dof;
cnt = round(t * freq) + 1;
if size(Td, 3) > 1
    Xd = Td(:,:, cnt);
else
    Xd = Td;
end
q = y(1:n);
qd = y(n + 1 : 2 * n);
q0 = inverse_kin_kuka_robot_kesai_near(robot, Xd, kesai, q);
tao = -Dn * qd - Kn * (q - q0');
[Jb, ~] = jacobian_matrix(robot, q);
M = mass_matrix(robot, q);
tao = tao - Jb' * ((Jb * (M \ Jb'))  \ (Jb * (M \ tao)));