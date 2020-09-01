function theta = R2EulXYZ(R)
% R = RxRyRz
% theta = [x,y,z]
r02 = R(1, 3);
r12 = R(2, 3);
r22 = R(3, 3);
r01 = R(1, 2);
r00 = R(1, 1);
r10 = R(2, 1);
r11 = R(2, 2);
theta = zeros(3, 1);
if r02 < 1
    if r02 > -1
        theta(1) = atan2(-r12, r22);
        theta(2) = asin(r02);
        theta(3) = atan2(-r01, r00);
    else
        theta(1) = -atan2(r10, r11);
        theta(2) = -pi / 2;
        theta(3) = 0;
    end
else
    theta(1) = atan2(r10, r11);
    theta(2) = pi / 2;
    theta(3) = 0;
end