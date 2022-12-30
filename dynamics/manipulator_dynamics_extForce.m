function yd = manipulator_dynamics_extForce(robot, tao, extForce, t, y)
% extForce is applied to the robot
% extForce(:,end) is applied to flange and others are applied to link frames
n = robot.dof;
yd = zeros(2*n,1);
q = y(1:n);
qd = y(n + 1 : end);
hqqd = gravity_velocity_torque(robot, q, qd);
yd(1:n) = qd;
J = jacobian_matrix_all(robot, q);
extf = extForce(t, y);
ext_torque = zeros(n,1);
for i = 1 : n
    if i < n
        Tbc = [eye(3), robot.com(i,:)'; 0 0 0 1];
        extf(:,i) = adjoint_T(Tbc)' * extf(:,i);
    end
    ext_torque = ext_torque + J(:,:,i)'*extf(:,i);
end
yd(n + 1 : end) = tao(t, y) - hqqd + ext_torque;