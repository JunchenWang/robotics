function [angles, flag] = inverse_kin_kuka_robot_kesai_near(robot, T, kesai, ref)
Td = T * tform_inv(robot.TCP);
tol = [1e-5, 1e-5];
[angles, flag] = inverse_kin_kuka_kesai_near(Td, kesai, ref, tol);