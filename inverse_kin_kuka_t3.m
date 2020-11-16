function angles = inverse_kin_kuka_t3(R, t, cfg, t3)
% inverse_kin_kuka kuka med的运动学逆解
% A3为冗余度信息, cfg符合kuka官方定义
% kesai redundancy
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

for i = [-1, 1]
    an = As(3,2); bn = Bs(3,2); cn = Cs(3,2);
    ad = -As(3,1); bd = -Bs(3,1); cd = -Cs(3,1);
    ap = i * ((cd-bd)*tan(t3)+(bn-cn));
    bp = 2*i*(ad*tan(t3)-an);
    cp = i*((bd+cd)*tan(t3)-(bn+cn));
    delta2 = bp^2 - 4*ap*cp;
    flag = 0;
    if delta2 < 0
        kesai = 0;
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
        break;
    else
        for kesai = [2*atan((-bp+sqrt(delta2)) / (2*ap)), 2*atan((-bp-sqrt(delta2)) / (2*ap))]
            t1 = atan2(i * (As(3,2) * sin(kesai) + Bs(3,2) * cos(kesai) + Cs(3,2)),...
                -i * (As(3,1) * sin(kesai) + Bs(3,1) * cos(kesai) + Cs(3,1)));
            if abs(t1-t3) < 1e-3
                theta2 = i * real(acos(As(3,3) * sin(kesai) + Bs(3,3) * cos(kesai) + Cs(3,3)));
                theta1_old = theta1;
                if abs(theta2) < eps2
                    % theta1 + theta3;
                    theta1 = atan2(As(2,1) * sin(kesai) + Bs(2,1) * cos(kesai) + Cs(2,1),...
                        As(1,1) * sin(kesai) + Bs(1,1) * cos(kesai) + Cs(1,1)) - t3;
                    theta3 = t3;
                else
                    theta1 = atan2(i * (As(2,3) * sin(kesai) + Bs(2,3) * cos(kesai) + Cs(2,3)),...
                        i * (As(1,3) * sin(kesai) + Bs(1,3) * cos(kesai) + Cs(1,3)));
                    theta3 = atan2(i * (As(3,2) * sin(kesai) + Bs(3,2) * cos(kesai) + Cs(3,2)),...
                        -i * (As(3,1) * sin(kesai) + Bs(3,1) * cos(kesai) + Cs(3,1)));
                end
                if cfg(1) == 1 && abs(theta1-theta1_old) <= pi / 2
                   flag = 1; break;
                elseif cfg(1) == -1 && abs(theta1-theta1_old) >= pi / 2
                   flag = 1; break;
                else
                    theta1 = theta1_old;
                end
            end
        end
        if flag == 1
            break;
        end
    end
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

