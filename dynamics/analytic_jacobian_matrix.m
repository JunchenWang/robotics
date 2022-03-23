function Ja = analytic_jacobian_matrix(Jb, T)
R = T(1:3,1:3);
r = logR(R)';
Ja = [inv(w_dr_A(r)), zeros(3); zeros(3), R] * Jb;