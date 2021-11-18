function gvtao = gravity_velocity_torque(mass, inertia, A, M, q, qd)
n = length(mass);
gvtao = inverse_dynamics(mass, inertia, A, M, eye(4), q, qd, zeros(n, 1), zeros(6, 1));

