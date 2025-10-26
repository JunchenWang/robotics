function [J2, J1] = derivative_r2r1(r2, r1)
q1 = r_q_converter(r1);
q2 = r_q_converter(r2);
q = ConcatenateQuaternions(q2, q1);
Jrq = derivative_r_q(q);
[Jqq2, Jqq1] = derivative_q2starq1(q2, q1);
Jq1r1 = derivative_q_r(r1);
Jq2r2 = derivative_q_r(r2);
J1 = Jrq * Jqq1 * Jq1r1;
J2 = Jrq * Jqq2 * Jq2r2;


