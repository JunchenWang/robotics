function q = EulerZYX2Quaternion(xyz)
xyz = xyz / 180 * pi;
cx = cos(xyz(1) / 2);
sx = sin(xyz(1) / 2);

cy = cos(xyz(2) / 2);
sy = sin(xyz(2) / 2);

cz = cos(xyz(3) / 2);
sz = sin(xyz(3) / 2);

q = [cx * cy * cz + sx * sy * sz
    sx * cy * cz - cx * sy * sz
    cx * sy * cz + sx * cy * sz
    cx * cy * sz - sx * sy * cz];