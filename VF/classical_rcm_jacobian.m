function [P2, rcm, Jrcm, J2] = classical_rcm_jacobian(robot, q, lambda, L)
Jb = jacobian_matrix(robot, q);
T = forward_kin_general(robot, q);
Tp1 = [eye(3), -[0, 0, 0]'; 0 0 0 1];
Tp2 = [eye(3), -[0, 0, L]'; 0 0 0 1];
R = T(1:3,1:3);
p = T(1:3,4);
J1 = adjoint_T(Tp1) * Jb;
J2 = adjoint_T(Tp2) * Jb;
J1 = R * J1(4:6,:);
J2 = R * J2(4:6,:);
P1 = p + R * [0, 0, 0]';
P2 = p + R * [0, 0, L]';
rcm = P1 + (P2 - P1) * lambda;
Jrcm = [J1 + lambda * (J2 - J1), P2 - P1];





