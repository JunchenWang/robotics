function [dJa, Ja, dT, T] = derivative_jacobian_matrix_analytic(robot, q, qd)
[dJb, Jb, dT, T] = derivative_jacobian_matrix(robot, q, qd);
if isrow(qd)
    qd = qd';
end
R = T(1:3,1:3);
r = logR(R);
A = w_dr_A(r);
Vb = Jb * qd;
wb = Vb(1:3);
dr =A \ wb;
dA = derivative_Ar(r, dr);
invA = inv(A);
Ja = [invA, zeros(3); zeros(3), R] * Jb;
dJa = [-invA * dA * invA, zeros(3); zeros(3), R * so_w(wb)] * Jb...
      + [invA, zeros(3); zeros(3), R] * dJb;