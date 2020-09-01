function q = RM2Quaternion(R)
r = AngleAxisFromRotation(R);
q = r_q_converter(r);