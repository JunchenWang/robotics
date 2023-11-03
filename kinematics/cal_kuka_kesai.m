function kesai = cal_kuka_kesai(angles)
d1 = .34;
d3 = .4;
d5 = .4;
d7 = .126;
if abs(angles(4)) < 1e-6 % 关节4为0，奇异位置
    kesai = 0;
    return;
elseif angles(4) > 0
    cfg4 = 1;
else
    cfg4 = -1;
end
eps1 = 1e-12;%not change
T = forward_kin_kuka(angles);
t = T(1:3,4);
R = T(1:3,1:3);
z = [0,0,1]';
T3 = forward_kin_kuka(angles(1:3));
E = T3(1:3,4);
S = [0, 0, d1]';
p02 = [0, 0, d1]';
p67 = [0, 0, d7]';
p26 = t - p02 - R * p67;
l_p26 = norm(p26);
p26_hat = p26 / l_p26;
l2_p26 = p26'*p26;
theta3 = 0;
if abs(abs(dot(p26_hat, z)) - 1) < eps1
    theta1 = 0;
else
    theta1 = atan2(p26(2), p26(1));
end
phi = real(acos((d3*d3 + l2_p26 - d5*d5)/(2*d3*l_p26)));
theta2 = atan2(sqrt(p26(1)^2 + p26(2)^2), p26(3)) + cfg4*phi;
T3_ = forward_kin_kuka([theta1, theta2, theta3]);
E_ = T3_(1:3,4);
SW = p26_hat;
SE_ = E_-S;
SE_ = SE_ / norm(SE_);
SE = E-S;
SE = SE / norm(SE);
% if abs(abs(dot(SE, SE_)) - 1) < eps1
%     kesai = 0;
%     return;
% end
x = cross(SW, SE_);
x = x / norm(x);
y = cross(SW, SE);
y = y / norm(y);
kesai = real(acos(x'*y));
if cross(x,y)'* SW < 0
    kesai = -kesai;
end
