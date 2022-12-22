function yd = manipulator_dynamics_extForce(robot, tao, extForce, t, y)
% extForce is applied to the robot
n = robot.dof;
yd = zeros(2*n,1);
q = y(1:n);
qd = y(n + 1 : end);
hqqd = gravity_velocity_torque(robot, q, qd);
yd(1:n) = qd;
J = jacobian_matrix_all(robot, q);
extf = extForce(t, y);
ext_torque = zeros(7,1);
for i = 1 : 7
    ext_torque = ext_torque + J(:,:,i)'*extf(:,i);
end
yd(n + 1 : end) = tao(t, y) - hqqd + ext_torque;