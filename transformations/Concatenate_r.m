function r = Concatenate_r(r2, r1)
% r = AngleAxisFromRotation(RotationByAxisAngleRep(r2) * RotationByAxisAngleRep(r1));
r = r_q_converter(ConcatenateQuaternions(r_q_converter(r2), r_q_converter(r1)));