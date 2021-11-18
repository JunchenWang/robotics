function vtao = velocity_torque(mass, inertia, A, M, q, qd)
gvtao = gravity_velocity_torque(mass, inertia, A, M, q, qd);
vtao = gvtao - gravity_torque(mass, inertia, A, M, q);
% M*qdd = b
% disp(Jb);
