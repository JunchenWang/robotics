function R = EulXYZ2R(theta)
% R = RxRyRz
cx = cos(theta(1));
sx = sin(theta(1));

cy = cos(theta(2));
sy = sin(theta(2));

cz = cos(theta(3));
sz = sin(theta(3));

R = [cy * cz, -cy * sz, sy;
    cz * sx * sy + cx * sz, cx * cz - sx * sy * sz, -cy * sx;
    -cx * cz * sy + sx * sz, cz * sx + cx * sy * sz, cx * cy];
