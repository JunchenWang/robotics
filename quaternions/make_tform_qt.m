function T = make_tform_qt(qt)
T = [RMByQuaternion(qt(1:4)), [qt(5); qt(6);qt(7)]; 0 0 0 1];