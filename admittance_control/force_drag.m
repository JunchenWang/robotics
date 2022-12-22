function speedj = force_drag(robot, Tsensor, Tcp, force_sensor, q, b, w)
Jb = jacobian_matrix(robot, q);
Jtcp = adjoint_T(InvertT(Tcp)) * Jb;
Ftcp = adjoint_T(InvertT(Tsensor) * Tcp)' * force_sensor .* w';
speedj = pinv(Jtcp) * Ftcp / b;
