function qdd = forward_dynamics_fext(robot, q, qd, tao, fext)
% fext is applied to the robot, different with forward_dynamics
n = robot.dof;
com = robot.com;
ME = robot.ME * robot.TCP;
[Mq, J] = mass_matrix(robot, q);
hqqd = inverse_dynamics(robot, q, qd, zeros(n, 1), zeros(6, 1));
ext_torque = zeros(n,1);
Tb = eye(4);
for i = 1 : n
    if i < n
       Tb(1:3,4) = com(i,:);
       J(:,:,i) = adjoint_T(Tb) * J(:,:,i);
    else
       J(:,:,n) = adjoint_T(tform_inv(ME)) * J(:,:,n);
    end
    ext_torque = ext_torque + J(:,:,i)'* fext(:,i);
end
b = tao - hqqd + ext_torque;
qdd = Mq \ b;
% M*qdd = b
% disp(Jb);
