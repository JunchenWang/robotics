function [q, flag] = inverse_kin_UR(Td, param, cfg, tol, ref6)
% param: |a2|, |a3|, d1, d4, d5, d6
% cfg: +- of q1, q3, q5
% tol: rotaton and position tolerance for last validation
% ref6: q6's angle if singularity
% 全部采用鲁棒性求解，例如real(acos(...))和atan2()，最后用正解验证
if nargin < 5
    ref6 = 0;
end
if nargin < 4
    tol = [1e-5, 1e-5];
end
if nargin < 3
    cfg = [1,1,1];
end
a2 = -param(1);
a3 = -param(2);
d1 = param(3);
d4 = param(4);
d5 = param(5);
d6 = param(6);

nx = Td(1,1);
ny = Td(2,1);
nz = Td(3,1);

ox = Td(1,2);
oy = Td(2,2);
oz = Td(3,2);

ax = Td(1,3);
ay = Td(2,3);
az = Td(3,3);

px = Td(1,4);
py = Td(2,4);
pz = Td(3,4);

q = zeros(6,1);
%% 计算q1
m1 = px - d6 * ax;
n1 = d6 * ay - py;
phi = atan2(m1, n1);
q(1) = mod(cfg(1) * real(acos(d4 / sqrt(m1^2 + n1^2))) + phi + pi, 2 * pi) - pi;
s1 = sin(q(1));
c1 = cos(q(1));
%% 计算q5
q(5) = cfg(3) * real(acos(s1 * ax - c1 * ay));
s5 = sin(q(5));
c5 = cos(q(5));

%% 计算q6
m6 = s1 * nx - c1 * ny;
n6 = c1 * oy - s1 * ox;
q(6) = atan2(n6 * s5, m6 * s5);
% 奇异值
if abs(s5) < 1e-7
    q(6) = ref6;
end
s6 = sin(q(6));
c6 = cos(q(6));
%% 计算q234
r31 = c5 * (c6 * nz - s6 * oz) - s5 * az;
r32 = s6 * nz + c6 * oz;
r14 = d5 * (s6 * (c1 * nx + s1 * ny) + c6 * (c1 * ox + s1 * oy)) - d6 * (c1 * ax + s1 * ay) + c1 * px + s1 * py;
r34 = d5 * (s6 * nz + c6 * oz) - d6 * az + pz - d1;
q234 = atan2(r31, r32);
%% 计算q3
q(3) = cfg(2) * real(acos((r14^2 + r34^2 - a2^2 - a3^2) / (2 * a2 * a3)));
s3 = sin(q(3));
c3 = cos(q(3));
%% 计算q2
m2 = a3 * s3;
n2 = a3 * c3 + a2;
q(2) = atan2(n2 * r34 - m2 * r14, m2 *r34 + n2 * r14);
%% 计算q4
q(4) = mod(q234 - q(2) - q(3) + pi, 2 * pi) - pi;
T = forward_kin_UR(q, param);
dR = norm(logR(T(1:3,1:3)' * Td(1:3,1:3)));
dp = norm(Td(1:3,4) - T(1:3,4));
if dR <= tol(1) && dp <= tol(2)
    flag = 1;
else
    flag = 0;
end
