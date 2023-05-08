function test_derivative_Ja
robot = convert_robot_tree(importrobot('E:\data\URDF\ur_5e-calibrated\ur_description\urdf\ur5e-A302.urdf'));
q = -rand(1,6);
qd = rand(1,6);
dJa = derivative_jacobian_matrix_analytic(robot, q, qd);
dt = 1e-6;
Ja1 = jacobian_matrix_analytic(robot, q);
Ja2 = jacobian_matrix_analytic(robot, q + dt * qd);
dJa2 = (Ja2 - Ja1) / dt;
norm(dJa)
norm(dJa2 - dJa)
