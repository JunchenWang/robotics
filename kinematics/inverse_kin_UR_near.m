function [q, flag] = inverse_kin_UR_near(Td, param, ref, tol)
q = ref(:);
flag = 0;
dist = 1e9;
for cfg1 = -1 : 2 : 1
    for cfg2 =  -1 : 2 : 1
        for cfg3 = -1 : 2 : 1
            cfg = [cfg1, cfg2, cfg3];
            [ang, f] = inverse_kin_UR(Td, param, cfg, tol, ref(6));
            if f == 1 && norm(ang-ref)<dist
                dist = norm(ang-ref);
                q = ang;
                flag = 1;
            end
        end
    end
end

% 以下代码和上面暴力便利的结果可能不同，局部最优不等于全局最优！
% if nargin < 4
%     tol = [1e-5, 1e-5];
% end
% 
% a2 = -param(1);
% a3 = -param(2);
% d1 = param(3);
% d4 = param(4);
% d5 = param(5);
% d6 = param(6);
% 
% nx = Td(1,1);
% ny = Td(2,1);
% nz = Td(3,1);
% 
% ox = Td(1,2);
% oy = Td(2,2);
% oz = Td(3,2);
% 
% ax = Td(1,3);
% ay = Td(2,3);
% az = Td(3,3);
% 
% px = Td(1,4);
% py = Td(2,4);
% pz = Td(3,4);
% 
% q = zeros(6,1);
% %% 计算q1
% m1 = px - d6 * ax;
% n1 = d6 * ay - py;
% phi = atan2(m1, n1);
% q1 = mod(real(acos(d4 / sqrt(m1^2 + n1^2))) + phi + pi, 2 * pi) - pi;
% q1_ = mod(-real(acos(d4 / sqrt(m1^2 + n1^2))) + phi + pi, 2 * pi) - pi;
% if abs(q1 - ref(1)) < abs(q1_ - ref(1))
%     q(1) = q1;
% else
%     q(1) = q1_;
% end
% s1 = sin(q(1));
% c1 = cos(q(1));
% %% 计算q5
% 
% q5 = real(acos(s1 * ax - c1 * ay));
% q5_ = -q5;
% if abs(q5 - ref(5)) < abs(q5_ - ref(5))
%     q(5) = q5;
% else
%     q(5) = q5_;
% end
% s5 = sin(q(5));
% c5 = cos(q(5));
% 
% %% 计算q6
% m6 = s1 * nx - c1 * ny;
% n6 = c1 * oy - s1 * ox;
% 
% if abs(s5) < 1e-7 % not change
%     q(6) = ref(6); %奇异位置
% else
%     q(6) = atan2(n6 * s5, m6 * s5);
% end
% s6 = sin(q(6));
% c6 = cos(q(6));
% %% 计算q234
% r31 = c5 * (c6 * nz - s6 * oz) - s5 * az;
% r32 = s6 * nz + c6 * oz;
% r14 = d5 * (s6 * (c1 * nx + s1 * ny) + c6 * (c1 * ox + s1 * oy)) - d6 * (c1 * ax + s1 * ay) + c1 * px + s1 * py;
% r34 = d5 * (s6 * nz + c6 * oz) - d6 * az + pz - d1;
% q234 = atan2(r31, r32);
% %% 计算q3
% q3 = real(acos((r14^2 + r34^2 - a2^2 - a3^2) / (2 * a2 * a3)));
% q3_ = -q3;
% if abs(q3 - ref(3)) < abs(q3_ - ref(3))
%     q(3) = q3;
% else
%     q(3) = q3_;
% end
% s3 = sin(q(3));
% c3 = cos(q(3));
% %% 计算q2
% m2 = a3 * s3;
% n2 = a3 * c3 + a2;
% q(2) = atan2(n2 * r34 - m2 * r14, m2 *r34 + n2 * r14);
% %% 计算q4
% q(4) = mod(q234 - q(2) - q(3) + pi, 2 * pi) - pi;
% T = forward_kin_UR(q, param);
% dR = norm(logR(T(1:3,1:3)' * Td(1:3,1:3)));
% dp = norm(Td(1:3,4) - T(1:3,4));
% if dR <= tol(1) && dp <= tol(2)
%     flag = 1;
% else
%     flag = 0;
% end