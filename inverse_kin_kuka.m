function angles = inverse_kin_kuka(R, t, cfg, kesai)
% inverse_kin_kuka kuka med的运动学逆解,冗余信息为kesai
% cfg signs of axis 2,4,6
% kesai redundancy
lower = [-170, -120, -170, -120 ,-170, -120, -175] / 180 * pi;
upper = -lower;
eps1 = 1e-6;
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
if norm(cross(p26_hat, z)) < eps1
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


kesai = anlysis_kesai(As, Bs, Cs, Aw, Bw, Cw, cfg, lower, upper, kesai);
kesai = kesai(1);
theta2 = cfg(1) * real(acos(As(3,3) * sin(kesai) + Bs(3,3) * cos(kesai) + Cs(3,3)));
if abs(theta2) < eps2
    % theta1 + theta3;
    theta1 = atan2(As(2,1) * sin(kesai) + Bs(2,1) * cos(kesai) + Cs(2,1),...
        As(1,1) * sin(kesai) + Bs(1,1) * cos(kesai) + Cs(1,1));
    theta3 = 0;
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
end


function [lb, ub] = theta2kesai(lower, upper, an, bn, cn, ad, bd, cd, cfg)
    at = cfg*(cn*bd-bn*cd);
    bt = cfg*(an*cd-cn*ad);
    ct = cfg*(an*bd-bn*ad);
    delta = at^2 + bt^2 - ct^2;
    if delta > 0
        kesai1 = 2*atan((at-sqrt(delta))/(bt-ct));
        kesai2 = 2*atan((at+sqrt(delta))/(bt-ct));
        theta1 = atan2(cfg*(an*sin(kesai1) + bn*cos(kesai1) + cn), cfg*(ad*sin(kesai1) + bd*cos(kesai1) + cd));
        theta2 = atan2(cfg*(an*sin(kesai2) + bn*cos(kesai2) + cn), cfg*(ad*sin(kesai2) + bd*cos(kesai2) + cd));
    else
        flag = 0;
    end
    theta = lower;
    ap = cfg * ((cd-bd)*tan(theta)+(bn-cn));
    bp = 2*cfg*(ad*tan(theta)-an);
    cp = cfg*((bd+cd)*tan(theta)-(bn+cn));
    delta = bp^2 - 4*ap*cp;
    kesais = [];
    cnt = 0;
    for kesai = [2*atan((-bp+sqrt(delta2)) / (2*ap)), 2*atan((-bp-sqrt(delta2)) / (2*ap))]
        t3_back = atan2(cfg(1) * (As(3,2) * sin(kesai) + Bs(3,2) * cos(kesai) + Cs(3,2)),...
            -cfg(1) * (As(3,1) * sin(kesai) + Bs(3,1) * cos(kesai) + Cs(3,1)));
        if abs(t3 - t3_back) < 1e-3
            cnt = cnt + 1;
            kesais(cnt) = kesai;
        end
    end
end

function kesais = anlysis_kesai(As, Bs, Cs, Aw, Bw, Cw, cfg, lower, upper, t3)
an = As(3,2); bn = Bs(3,2); cn = Cs(3,2);
ad = -As(3,1); bd = -Bs(3,1); cd = -Cs(3,1);
at = cfg(1)*(cn*bd-bn*cd);
bt = cfg(1)*(an*cd-cn*ad);
ct = cfg(1)*(an*bd-bn*ad);
delta = at^2 + bt^2 - ct^2;

kesais = [0, 0];
y = @(kesai) atan2(cfg(1) * (As(3,2) .* sin(kesai) + Bs(3,2) .* cos(kesai) + Cs(3,2)),...
    -cfg(1) * (As(3,1) .* sin(kesai) + Bs(3,1) .* cos(kesai) + Cs(3,1)));
x = linspace(-pi, pi, 100);
plot(x, y(x));
theta = t3;
 ap = cfg(1) * ((cd-bd)*tan(theta)+(bn-cn));
    bp = 2*cfg(1)*(ad*tan(theta)-an);
    cp = cfg(1)*((bd+cd)*tan(theta)-(bn+cn));
    delta2 = bp^2 - 4*ap*cp;
    kesai = [2*atan((-bp+sqrt(delta2)) / (2*ap)), 2*atan((-bp-sqrt(delta2)) / (2*ap))]
end
