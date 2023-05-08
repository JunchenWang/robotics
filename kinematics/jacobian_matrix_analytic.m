function [Ja, T] = jacobian_matrix_analytic(robot, q)
% (dr,dp)
[Jb, T] = jacobian_matrix(robot, q);
R = T(1:3,1:3);
r = logR(R);
A = w_dr_A(r);
Ja = [inv(A), zeros(3); zeros(3), R] * Jb;