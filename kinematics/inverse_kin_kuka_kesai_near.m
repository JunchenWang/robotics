function [angles, flag] = inverse_kin_kuka_kesai_near(Td, kesai, ref, tol)
angles = ref(:);
dist = 1e12;
flag = 0;
for cfg1 = -1 : 2 : 1
    for cfg2 =  -1 : 2 : 1
        for cfg3 = -1 : 2 : 1
            cfg = [cfg1, cfg2, cfg3];
            [ang, f] = inverse_kin_kuka_kesai(Td, cfg, kesai, tol, ref);
            if f && norm(ang-ref)<dist
                dist = norm(ang-ref);
                angles = ang;
                flag = 1;
            end
        end
    end
end

% if nargin < 4
%     tol = [1e-5, 1e-5];
% end
% cfg = [1,1,1];
% eps1 = 1e-12;%not change
% eps0 = 1e-8;%not change
% z = [0, 0, 1]';
% d1 = .34;
% d3 = .4;
% d5 = .4;
% d7 = .126;
% R = Td(1:3,1:3);
% t = Td(1:3, 4);
% p02 = [0, 0, d1]';
% p67 = [0, 0, d7]';
% p26 = t - p02 - R * p67;
% l_p26 = norm(p26);
% p26_hat = p26 / l_p26;
% l2_p26 = p26'*p26;
% theta3 = 0;
% L = (l2_p26 - d3*d3 - d5*d5) / (2*d3*d5);
% if ref(4) < 0
%     cfg(2) = -1;
% end
% theta4 = cfg(2) * real(acos(L));
% % R34 make both shoulder and wrist a zyz spherical joint
% R34 = [cos(theta4), 0, -sin(theta4);0, 1, 0;sin(theta4), 0, cos(theta4)];
% % theta1 = atan2(p26(2), p26(1));
% if abs(abs(dot(p26_hat, z)) - 1) < eps1
%     theta1 = 0;
% else
%     theta1 = atan2(p26(2), p26(1));
% end
% phi = real(acos((d3*d3 + l2_p26 - d5*d5)/(2*d3*l_p26)));
% theta2 = atan2(sqrt(p26(1)^2 + p26(2)^2), p26(3)) + cfg(2)*phi;
% T03 = forward_kin_kuka([theta1, theta2, theta3]);
% R03 = T03(1:3, 1:3);
% % shoulder
% As = so_w(p26_hat) * R03;
% Bs = -so_w(p26_hat) * As;
% Cs = p26_hat * p26_hat' * R03;
% % wrist
% Aw = R34' * As' * R;
% Bw = R34' * Bs' * R;
% Cw = R34' * Cs' * R;
% if ref(2) < 0
%     cfg(1) = -1;
% end
% t2 = As(3,3) * sin(kesai) + Bs(3,3) * cos(kesai) + Cs(3,3);
% theta2 = cfg(1) * real(acos(t2));
% if abs(abs(t2)-1) < eps0
%     theta1 = ref(1);
%     if t2 > 0
%         % theta1 + theta3;
%         theta1and3 = atan2(As(2,1) * sin(kesai) + Bs(2,1) * cos(kesai) + Cs(2,1),...
%             As(1,1) * sin(kesai) + Bs(1,1) * cos(kesai) + Cs(1,1));
%         theta3= mod(theta1and3-theta1 + pi, 2 * pi) - pi;
%     else
%         %theta3 - theta1;
%         theta3_1 = atan2(As(1,2) * sin(kesai) + Bs(1,2) * cos(kesai) + Cs(1,2),...
%             As(2,2) * sin(kesai) + Bs(2,2) * cos(kesai) + Cs(2,2));
%         theta3 = mod(theta3_1 + theta1 + pi, 2 * pi) - pi;
%     end
% else
%     theta1 = atan2(cfg(1) * (As(2,3) * sin(kesai) + Bs(2,3) * cos(kesai) + Cs(2,3)),...
%         cfg(1) * (As(1,3) * sin(kesai) + Bs(1,3) * cos(kesai) + Cs(1,3)));
%     theta3 = atan2(cfg(1) * (As(3,2) * sin(kesai) + Bs(3,2) * cos(kesai) + Cs(3,2)),...
%         -cfg(1) * (As(3,1) * sin(kesai) + Bs(3,1) * cos(kesai) + Cs(3,1)));
% end
% t6 = Aw(3,3) * sin(kesai) + Bw(3,3) * cos(kesai) + Cw(3,3);
% if ref(6) < 0
%     cfg(3) = -1;
% end
% theta6 = cfg(3) * real(acos(t6));
% if abs(abs(t6) - 1) < eps0
%     theta5 = ref(5);
%     if t6 > 0
%         % theta5 + theta7 ;
%         theta5and7 = atan2(Aw(2,1) * sin(kesai) + Bw(2,1) * cos(kesai) + Cw(2,1),...
%             Aw(1,1) * sin(kesai) + Bw(1,1) * cos(kesai) + Cw(1,1));
%         theta7= mod(theta5and7-theta5 + pi, 2 * pi) - pi;
%     else
%         % theta7 - theta5 ;
%         theta7_5 = atan2(Aw(1,2) * sin(kesai) + Bw(1,2) * cos(kesai) + Cw(1,2),...
%             Aw(2,2) * sin(kesai) + Bw(2,2) * cos(kesai) + Cw(2,2));
%         theta7 = mod(theta7_5 + theta5 + pi, 2 * pi) - pi;
%     end
% else
%     theta5 = atan2(cfg(3) * (Aw(2,3) * sin(kesai) + Bw(2,3) * cos(kesai) + Cw(2,3)),...
%         cfg(3) * (Aw(1,3) * sin(kesai) + Bw(1,3) * cos(kesai) + Cw(1,3)));
%     theta7 = atan2(cfg(3) * (Aw(3,2) * sin(kesai) + Bw(3,2) * cos(kesai) + Cw(3,2)),...
%         -cfg(3) * (Aw(3,1) * sin(kesai) + Bw(3,1) * cos(kesai) + Cw(3,1)));
% end
% angles = [theta1, theta2, theta3, theta4, theta5, theta6, theta7]';
% T = forward_kin_kuka(angles);
% dR = norm(logR(T(1:3,1:3)' * Td(1:3,1:3)));
% dp = norm(Td(1:3,4) - T(1:3,4));
% if dR <= tol(1) && dp <= tol(2)
%     flag = 1;
% else
%     flag = 0;
% end
% end
