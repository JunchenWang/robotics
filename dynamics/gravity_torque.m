function gtao = gravity_torque(mass, inertia, A, M, q)
n = length(mass);
gtao = inverse_dynamics(mass, inertia, A, M, eye(4), q, zeros(n, 1), zeros(n, 1), zeros(6, 1));
