function angles = inverse_kin_kuka(R, t, cfg, kesai)
% inverse_kin_kuka kuka med的运动学逆解
% cfg signs of axis 2,4,6
% kesai redundancy
eps1 = 1e-9;
eps2 = 1e-6;
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
theta4 = cfg(2) * acos((l2_p26 - d3*d3 - d5*d5) / (2*d3*d5));
if abs(cross(p26, z)) < eps1
    theta1 = 0;
else
    theta1 = atan2(p26(2), p26(1));
end
phi = acos((d3*d3 + l2_p26 - d5*d5)/(2*d3*l_p26));
theta2 = atan2(sqrt(p26(1)^2 + p26(2)^2), p26(3)) + cfg(2)*phi;
T03 = forward_kin_kuka([theta1, theta2, theta3]);
R03 = T03(1:3, 1:3);
A = so_w(p26_hat) * R03;
B = -so_w(p26_hat) * A;
C = p26_hat * p26_hat' * R03;
theta2 = cfg(1) * real(acos(A(3,3) * sin(kesai) + B(3,3) * cos(kesai) + C(3,3)));
if abs(theta2) < eps2 
    % theta1 + theta3;
    theta1 = atan2(A(2,1) * sin(kesai) + B(2,1) * cos(kesai) + C(2,1),...
                 A(1,1) * sin(kesai) + B(1,1) * cos(kesai) + C(1,1));
    theta3 = 0;
else
theta1 = atan2(cfg(1) * (A(2,3) * sin(kesai) + B(2,3) * cos(kesai) + C(2,3)),...
               cfg(1) * (A(1,3) * sin(kesai) + B(1,3) * cos(kesai) + C(1,3)));
theta3 = atan2(cfg(1) * (A(3,2) * sin(kesai) + B(3,2) * cos(kesai) + C(3,2)),...
               -cfg(1) * (A(3,1) * sin(kesai) + B(3,1) * cos(kesai) + C(3,1)));
an = A(3,2); bn = B(3,2); cn = C(3,2);
ad = -A(3,1); bd = -B(3,1); cd = -C(3,1);
at = cfg(1)*(cn*bd-bn*cd);
bt = cfg(1)*(an*cd-cn*ad);
ct = cfg(1)*(an*bd-bn*ad);
delta = at^2 + bt^2 - ct^2;
y = @(kesai) atan2(cfg(1) * (A(3,2) .* sin(kesai) + B(3,2) .* cos(kesai) + C(3,2)),...
               -cfg(1) * (A(3,1) .* sin(kesai) + B(3,1) .* cos(kesai) + C(3,1)));
% x = linspace(-pi, pi, 100);
% plot(x, y(x));
ap = cfg(1) * ((cd-bd)*tan(theta3)+(bn-cn));
bp = 2*cfg(1)*(ad*tan(theta3)-an);
cp = cfg(1) * ((bd+cd)*tan(theta3)-(bn+cn));
delta2 = bp^2 - 4*ap*cp;
disp(theta3 / pi * 180);
kesai2 = 2*atan((-bp+sqrt(delta2)) / (2*ap)) / pi * 180;
kesai2_1 = 2*atan((-bp-sqrt(delta2)) / (2*ap)) / pi * 180;
end
R34 = [cos(theta4), 0, -sin(theta4);
       0, 1, 0;
       sin(theta4), 0, cos(theta4)];
A = R34' * A' * R;
B = R34' * B' * R;
C = R34' * C' * R;
theta6 = cfg(3) * real(acos(A(3,3) * sin(kesai) + B(3,3) * cos(kesai) + C(3,3)));
if abs(theta6) < eps2
    % theta5 + theta7 ;
    theta7 = 0;
    theta5 = atan2(A(2,1) * sin(kesai) + B(2,1) * cos(kesai) + C(2,1),...
                 A(1,1) * sin(kesai) + B(1,1) * cos(kesai) + C(1,1));
else
theta5 = atan2(cfg(3) * (A(2,3) * sin(kesai) + B(2,3) * cos(kesai) + C(2,3)),...
               cfg(3) * (A(1,3) * sin(kesai) + B(1,3) * cos(kesai) + C(1,3)));
theta7 = atan2(cfg(3) * (A(3,2) * sin(kesai) + B(3,2) * cos(kesai) + C(3,2)),...
               -cfg(3) * (A(3,1) * sin(kesai) + B(3,1) * cos(kesai) + C(3,1)));
end
angles = [theta1, theta2, theta3, theta4, theta5, theta6, theta7] / pi * 180;

% T04 = forward_kin_kuka([theta1, theta2, theta3, theta4]);
% R04 = T04(1:3,1:3);
% p04 = T04(1:3, 4);
% p06 = p02 + p26;
% p24 = p04 - p02;
% v_sew = cross(p24/norm(p24), p26 / l_p26);
% Rk = exp_w(p26 / l_p26 * kesai);

end

