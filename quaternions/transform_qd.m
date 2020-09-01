function [R t] = transform_qd(qd)
q0 = qd(1 : 4);
R = RMByQuaternion(q0);
t = 2 * ConcatenateQuaternions(qd(5 : 8), Conjugate_q(q0));
t = t(2 : end);