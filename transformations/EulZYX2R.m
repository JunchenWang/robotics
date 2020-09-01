function R = EulZYX2R(theta)
% theta: in the order of z,y,x
cy = cos(theta(2));
cz = cos(theta(1));
cx = cos(theta(3));

sy = sin(theta(2));
sz = sin(theta(1));
sx = sin(theta(3));

R = [cy * cz, cz * sx * sy - cx * sz, cx * cz * sy + sx * sz;
     cy * sz, cx * cz + sx * sy * sz, -cz * sx + cx * sy * sz;
     -sy, cy * sx, cx * cy];
