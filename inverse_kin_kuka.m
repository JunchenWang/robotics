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
R34 = [cos(theta4), 0, -sin(theta4);0, 1, 0;sin(theta4), 0, cos(theta4)];
if abs(cross(p26, z)) < eps1
    theta1 = 0;
else
    theta1 = atan2(p26(2), p26(1));
end
phi = acos((d3*d3 + l2_p26 - d5*d5)/(2*d3*l_p26));
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

t3 = kesai;
% kesai = analyze_kesai(cfg, As, Bs, Cs, Aw, Bw, Cw, t3);
kesai = theta2kesai(cfg(1), As, Bs, Cs, t3);
theta2 = cfg(1) * real(acos(As(3,3) * sin(kesai) + Bs(3,3) * cos(kesai) + Cs(3,3)));
if abs(theta2) < eps2
    % theta1 + theta3;
    theta1 = atan2(As(2,1) * sin(kesai) + Bs(2,1) * cos(kesai) + Cs(2,1),...
        As(1,1) * sin(kesai) + Bs(1,1) * cos(kesai) + Cs(1,1)) - t3;
    theta3 = t3;
else
    theta1 = atan2(cfg(1) * (As(2,3) * sin(kesai) + Bs(2,3) * cos(kesai) + Cs(2,3)),...
        cfg(1) * (As(1,3) * sin(kesai) + Bs(1,3) * cos(kesai) + Cs(1,3)));
    theta3 = atan2(cfg(1) * (As(3,2) * sin(kesai) + Bs(3,2) * cos(kesai) + Cs(3,2)),...
        -cfg(1) * (As(3,1) * sin(kesai) + Bs(3,1) * cos(kesai) + Cs(3,1)));
end
theta6 = cfg(3) * real(acos(Aw(3,3) * sin(kesai) + Bw(3,3) * cos(kesai) + Cw(3,3)));
if abs(theta6) < eps2
    % theta5 + theta7 ;
    theta7 = 0;
    theta5 = atan2(Aw(2,1) * sin(kesai) + Bw(2,1) * cos(kesai) + Cw(2,1),...
        Aw(1,1) * sin(kesai) + Bw(1,1) * cos(kesai) + Cw(1,1));
else
    theta5 = atan2(cfg(3) * (Aw(2,3) * sin(kesai) + Bw(2,3) * cos(kesai) + Cw(2,3)),...
        cfg(3) * (Aw(1,3) * sin(kesai) + Bw(1,3) * cos(kesai) + Cw(1,3)));
    theta7 = atan2(cfg(3) * (Aw(3,2) * sin(kesai) + Bw(3,2) * cos(kesai) + Cw(3,2)),...
        -cfg(3) * (Aw(3,1) * sin(kesai) + Bw(3,1) * cos(kesai) + Cw(3,1)));
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

function kesai = theta2kesai(cfg, A, B, C, t)
an = A(3,2); bn = B(3,2); cn = C(3,2);
ad = -A(3,1); bd = -B(3,1); cd = -C(3,1);
at = cfg*(cn*bd-bn*cd);
bt = cfg*(an*cd-cn*ad);
ct = cfg*(an*bd-bn*ad);
delta = at^2 + bt^2 - ct^2;
y = @(kesai) atan2(cfg(1) * (A(3,2) .* sin(kesai) + B(3,2) .* cos(kesai) + C(3,2)),...
    -cfg(1) * (A(3,1) .* sin(kesai) + B(3,1) .* cos(kesai) + C(3,1)));
x = linspace(-pi, pi, 100);
plot(x, y(x));
ap = cfg * ((cd-bd)*tan(t)+(bn-cn));
bp = 2*cfg*(ad*tan(t)-an);
cp = cfg*((bd+cd)*tan(t)-(bn+cn));
delta2 = bp^2 - 4*ap*cp;
% if abs(ap) < 1e-9
%     kesai = 2*atan(-cp/bp);
if delta2 < 0
    kesai = 0;
else
    kesai_1 = 2*atan((-bp+sqrt(delta2)) / (2*ap));
    kesai_2 = 2*atan((-bp-sqrt(delta2)) / (2*ap));
    
    t1 = atan2(cfg * (A(3,2) * sin(kesai_1) + B(3,2) * cos(kesai_1) + C(3,2)),...
        -cfg * (A(3,1) * sin(kesai_1) + B(3,1) * cos(kesai_1) + C(3,1)));
    t2 = atan2(cfg * (A(3,2) * sin(kesai_2) + B(3,2) * cos(kesai_2) + C(3,2)),...
        -cfg * (A(3,1) * sin(kesai_2) + B(3,1) * cos(kesai_2) + C(3,1)));
    if abs(t-t1) < 1e-3
        kesai = kesai_1;
    elseif abs(t-t2) < 1e-3
        kesai = kesai_2;
    else
        kesai = 0;
    end
end
end

function kesai = analyze_kesai(cfg, As, Bs, Cs, Aw, Bw, Cw, theta3)
an = As(3,2); bn = Bs(3,2); cn = Cs(3,2);
ad = -As(3,1); bd = -Bs(3,1); cd = -Cs(3,1);
at = cfg(1)*(cn*bd-bn*cd);
bt = cfg(1)*(an*cd-cn*ad);
ct = cfg(1)*(an*bd-bn*ad);
delta = at^2 + bt^2 - ct^2;
y = @(kesai) atan2(cfg(1) * (As(3,2) .* sin(kesai) + Bs(3,2) .* cos(kesai) + Cs(3,2)),...
    -cfg(1) * (As(3,1) .* sin(kesai) + Bs(3,1) .* cos(kesai) + Cs(3,1)));
x = linspace(-pi, pi, 100);
plot(x, y(x));
ap = cfg(1) * ((cd-bd)*tan(theta3)+(bn-cn));
bp = 2*cfg(1)*(ad*tan(theta3)-an);
cp = cfg(1) * ((bd+cd)*tan(theta3)-(bn+cn));
delta2 = bp^2 - 4*ap*cp;
% if abs(ap) < 1e-9
%     kesai = 2*atan(-cp/bp);
if delta2 < 0
    kesai = 0;
else
    kesai_1 = 2*atan((-bp+sqrt(delta2)) / (2*ap));
    kesai_2 = 2*atan((-bp-sqrt(delta2)) / (2*ap));
    
    theta3_back = atan2(cfg(1) * (As(3,2) * sin(kesai_1) + Bs(3,2) * cos(kesai_1) + Cs(3,2)),...
        -cfg(1) * (As(3,1) * sin(kesai_1) + Bs(3,1) * cos(kesai_1) + Cs(3,1)));
    if abs(theta3-theta3_back) < 1e-3
        kesai = kesai_1;
    else
        kesai = kesai_2;
    end
end
end

