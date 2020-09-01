function theta = R2EulZYX(R)
% theta = [Z ,Y, X]
r20 = R(3, 1);
r10 = R(2, 1);
r00 = R(1, 1);
r21 = R(3, 2);
r22 = R(3, 3);
r12 = R(2, 3);
r11 = R(2, 2);

if r20 < 1
    if r20 > -1
        theta = [atan2(r10, r00); asin(-r20); atan2(r21, r22)];
    else
        theta = [-atan2(-r12, r11); pi / 2; 0];
    end
else
    theta = [atan2(-r12, r11); -pi / 2; 0];
end