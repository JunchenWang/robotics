function xyz = Quaternion2EulerZYX(q)

test = 2 * (q(1) * q(3) - q(4) * q(2));
%Gimbal lock check
if test > 1 - 1e-10
    xyz(1) = 0;
    xyz(2) = pi / 2;
    xyz(3) = -2 * atan2(q(2), q(1));
elseif test < 1e-10 - 1
    xyz(1) = 0;
    xyz(2) = -pi / 2;
    xyz(3) = 2 * atan2(q(2), q(1));
else
    xyz(1) = atan2(2 * (q(1) * q(2) + q(3) * q(4)), 1 - 2 * (q(2)^2 + q(3)^2));
    xyz(2) = asin(test);
    xyz(3) = atan2(2 * (q(1) * q(4) + q(2) * q(3)), 1 - 2 * (q(3)^2 + q(4)^2));
end
xyz = xyz' / pi * 180;
% R = RMByQuaternion(q);
% xyz = extractEulerZYX(R);