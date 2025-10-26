function [Jx Jb] = derivative_qdxqdbqdx_(qd_x, qd_b)
q_x = qd_x(1:4);
qe_x = qd_x(5:8);
q_b = qd_b(1:4);
qe_b = qd_b(5:8);

[Jqx1 Jqb1] = derivative_qxqbqx_(q_x, q_b);
[Jqx2 Jqb2] = derivative_qxqbqx_(q_x, qe_b);
[J1 J2] = derivative_q2starq1(q_x, ConcatenateQuaternions(q_b, Conjugate_q(qe_x)));
[J3 J4] = derivative_q2starq1(ConcatenateQuaternions(qe_x, q_b), Conjugate_q(q_x));

[J5 J6] = derivative_q2starq1(ConcatenateQuaternions(q_x, q_b), Conjugate_q(qe_x));
[J7 J8] = derivative_q2starq1(qe_x, ConcatenateQuaternions(q_b, Conjugate_q(q_x)));
I = -eye(4);
I(1,1) = 1;
Jx = [Jqx1, zeros(4,4);J1 + Jqx2 + J4 * I, J6 * I + J7];
J9 = derivative_q2starq1(q_b, Conjugate_q(qe_x));
J10 = derivative_q2starq1(q_b, Conjugate_q(q_x));
Jb = [Jqb1, zeros(4,4);  J2 * J9 + J8 * J10, Jqb2];
