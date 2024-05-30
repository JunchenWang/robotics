function test_orientation_jacobian


robot = convert_robot_tree2(importrobot('urdf\iiwa7\iiwa7.urdf'));

q = rand(7,1);
dq = rand(7,1);

[J, dJ, r] = orientation_jacobian(robot, q, dq);

dt = 1e-6;
q1 = q + dq * dt;
[J1, ~, r1] = orientation_jacobian(robot, q1, dq);

dr_truth = (r1 - r) / dt;
dJ_truth = (J1 - J) / dt;
disp(norm(dr_truth - J * dq));
disp(norm(dJ_truth - dJ));