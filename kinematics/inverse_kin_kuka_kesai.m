function angles = inverse_kin_kuka_kesai(R, t, cfg, kesai, hint)
% inverse_kin_kuka kuka med的运动学逆解,冗余信息为kesai
% cfg signs of axis 2,4,6
% kesai redundancy
% hint last angles
if nargin < 5
    hint = [];
end
angles = [];
eps1 = 1e-12;%not change
eps0 = 1e-7;%not change
z = [0, 0, 1]';
d1 = 340;
d3 = 400;
d5 = 400;
d7 =126;

p02 = [0, 0, d1]';
p67 = [0, 0, d7]';
p26 = t - p02 - R * p67;
l_p26 = norm(p26);
p26_hat = p26 / l_p26;
l2_p26 = p26'*p26;
theta3 = 0;
L = (l2_p26 - d3*d3 - d5*d5) / (2*d3*d5);
if abs(L) > 1
    return;
end
theta4 = cfg(2) * real(acos(L));
% R34 make both shoulder and wrist a zyz spherical joint
R34 = [cos(theta4), 0, -sin(theta4);0, 1, 0;sin(theta4), 0, cos(theta4)];
% theta1 = atan2(p26(2), p26(1));
if abs(abs(dot(p26_hat, z)) - 1) < eps1
    theta1 = 0;
else
    theta1 = atan2(p26(2), p26(1));
end
phi = real(acos((d3*d3 + l2_p26 - d5*d5)/(2*d3*l_p26)));
theta2 = atan2(sqrt(p26(1)^2 + p26(2)^2), p26(3)) + cfg(2)*phi;
T03 = forward_kin_kuka([theta1, theta2, theta3]);
R03 = T03(1:3, 1:3);
% shoulder
As = so_w(p26_hat) * R03;
Bs = -so_w(p26_hat) * As;
Cs = p26_hat * p26_hat' * R03;
% wrist
Aw = R34' * As' * R;
Bw = R34' * Bs' * R;
Cw = R34' * Cs' * R;

t2 = As(3,3) * sin(kesai) + Bs(3,3) * cos(kesai) + Cs(3,3);
theta2 = cfg(1) * real(acos(t2));
if abs(abs(t2)-1) < eps0
    if isempty(hint)
        theta1 = 0;
    else
        theta1 = hint(1);
    end
    if t2 > 0
        % theta1 + theta3;
        theta1and3 = atan2(As(2,1) * sin(kesai) + Bs(2,1) * cos(kesai) + Cs(2,1),...
            As(1,1) * sin(kesai) + Bs(1,1) * cos(kesai) + Cs(1,1));
        theta3=theta1and3-theta1;
    else
        %theta3 - theta1;
        theta3_1 = atan2(As(1,2) * sin(kesai) + Bs(1,2) * cos(kesai) + Cs(1,2),...
            As(2,2) * sin(kesai) + Bs(2,2) * cos(kesai) + Cs(2,2));
        theta3 = theta3_1 + theta1;
    end
else
    theta1 = atan2(cfg(1) * (As(2,3) * sin(kesai) + Bs(2,3) * cos(kesai) + Cs(2,3)),...
        cfg(1) * (As(1,3) * sin(kesai) + Bs(1,3) * cos(kesai) + Cs(1,3)));
    theta3 = atan2(cfg(1) * (As(3,2) * sin(kesai) + Bs(3,2) * cos(kesai) + Cs(3,2)),...
        -cfg(1) * (As(3,1) * sin(kesai) + Bs(3,1) * cos(kesai) + Cs(3,1)));
end
t6 = Aw(3,3) * sin(kesai) + Bw(3,3) * cos(kesai) + Cw(3,3);
theta6 = cfg(3) * real(acos(t6));
if abs(abs(t6) - 1) < eps0
    if isempty(hint)
        theta5 = 0;
    else
        theta5 = hint(5);
    end
    if t6 > 0
        % theta5 + theta7 ;
        theta5and7 = atan2(Aw(2,1) * sin(kesai) + Bw(2,1) * cos(kesai) + Cw(2,1),...
            Aw(1,1) * sin(kesai) + Bw(1,1) * cos(kesai) + Cw(1,1));
        theta7=theta5and7-theta5;
    else
        % theta7 - theta5 ;
        theta7_5 = atan2(Aw(1,2) * sin(kesai) + Bw(1,2) * cos(kesai) + Cw(1,2),...
            Aw(2,2) * sin(kesai) + Bw(2,2) * cos(kesai) + Cw(2,2));
        theta7 = theta7_5 + theta5;
    end
else
    theta5 = atan2(cfg(3) * (Aw(2,3) * sin(kesai) + Bw(2,3) * cos(kesai) + Cw(2,3)),...
        cfg(3) * (Aw(1,3) * sin(kesai) + Bw(1,3) * cos(kesai) + Cw(1,3)));
    theta7 = atan2(cfg(3) * (Aw(3,2) * sin(kesai) + Bw(3,2) * cos(kesai) + Cw(3,2)),...
        -cfg(3) * (Aw(3,1) * sin(kesai) + Bw(3,1) * cos(kesai) + Cw(3,1)));
end
angles = [theta1, theta2, theta3, theta4, theta5, theta6, theta7];
end

