function ext_torque = get_ext_torque(robot, q, fext)
% Fext 6-by-n, each column is external wrench imposed on the link of the robot 
J = jacobian_matrix_all(robot, q);
n = robot.dof;
ext_torque = zeros(n,1);
for i = 1 : n
    ext_torque = ext_torque + J(:,:,i)'* fext(:,i);
end