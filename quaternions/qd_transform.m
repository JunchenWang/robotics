function qd = qd_transform(R, t)
q0 = RM2Quaternion(R);
t = t(:);
qe = ConcatenateQuaternions([0; t / 2], q0);
qd = [ q0; qe];