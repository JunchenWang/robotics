function J = derivative_f_qd(qd)
q = qd(1:4);
qe = qd(5:8);
J = zeros(6, 8);

Jq = derivative_r_q(q);
J(1:3, :) = [Jq, zeros(3, 4)];
[J1 J2] = derivative_q2starq1(qe, Conjugate_q(q));
I = -eye(4);
I(1,1) = 1;
J3 = 2 * [J2 * I, J1];
J(4:6, :) = J3(2:4,:);