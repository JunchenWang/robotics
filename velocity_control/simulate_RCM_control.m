function simulate_RCM_control
% 双环控制7轴机器人轨迹，状态空间4个变量：位置，速度，速度误差积分，位置误差积分
port = udpport("byte");
opts = odeset('OutputFcn', @(t, y, flag) odeplot(t, y, flag, port));
robot = convert_robot_tree(importrobot('urdf\iiwa7\iiwa7.urdf'));
n = robot.dof;
y0 = zeros(4 * n, 1);

Kp_s = [100,100,100,100,100,100,100]';
Ki_s = [50,50,50,50,100,100,100]';

Kp_p = 20 * [1, 1, 1, 1, 1, 1]';
Ki_p = 0 * [1,1,1,1,1,1]';
Kd_p = 0 * [1,1,1,1,1,1]';

d = @(t, y) [7;6;5;4;3;2;1]*10;
J = [8;5;4;3;2;2;2];
B = 2;
r = 200;
tspan = [0, 10];

y0(1:n) = [0 75 0 -94 0 -81 0] / 180 * pi;
% y0(1:n) = [0 40 0 -80  -10 45 0] / 180 * pi;
ptp(port, y0(1:n)');
kesai = cal_kuka_kesai(y0);
Ts = forward_kin_general(robot, y0);
Rs = Ts(1:3,1:3);
ps = Ts(1:3,4);
% RCM param
p1 = 0;
p2 = 0.3;

lambda0 = 0.5;

P1 = ps + Rs * [0, 0, p1]';
P2 = ps + Rs * [0, 0, p2]';
Prcm = P1 + (P2 - P1) * lambda0;

y0(end) = lambda0;

p = @(t) desired_rcm_pos(t, Ts, p1, p2, lambda0);
v = @(t, y) task_rcm_controller(t, y, p, robot, Ts, p1, p2, kesai, Kp_p, Ki_p, Kd_p);
u = @(t, y) speed_controller(t, y, v, Kp_s, Ki_s);
dynamic = @(t, y) joint_motor_dynamic(t, y, u, d, J, B, r);

[t, y] = ode45(dynamic, tspan, y0, opts);

pos_d = zeros(3, length(t));
pos_actual = zeros(3, length(t));
P1 = zeros(3, length(t));
P2 = zeros(3, length(t));
dist = zeros(length(t),1);
rcm_error = zeros(2,length(t));
lambdas = zeros(length(t),1);
for i = 1 : length(t)
    [xyz, prcm] = p(t(i));
    pos_d(:,i) = ps + Rs * xyz;
    T = forward_kin_general(robot, y(i,1:7));
    lambdas(i) = y(i,end);
    R = T(1:3,1:3);
    pp = T(1:3,4);
    P1(:,i) = pp + R * [0, 0, p1]';
    P2(:,i) = pp + R * [0, 0, p2]';
    pos_actual(:,i) = P2(:,i);
    dist(i) = dist2line(prcm, P1(:,i), P2(:,i));
    rcm_error(1,i) = norm(y(i, 3*n + 4:3*n + 6));
    rcm_error(2,i) = norm(y(i, 3*n + 1:3*n + 3));
end
figure;
plot(t, dist,'LineWidth', 2);
% plot(t, dist, t, rcm_error(1,:), '--', t, rcm_error(2,:), '-*','LineWidth', 2);
xlabel("$t$/s", 'interpreter','latex');
ylabel('RCM误差/m', 'interpreter','latex');
% yticks([0,.2, .4, .6, .8, 1.0, 1.2, 1.4]);
set(gca,'FontSize', 36);
% figure;
% plot(t, pos_actual(1:3,:), t, pos_d(1:3,:), '--', 'LineWidth', 2);
% xlabel("$t$/s", 'interpreter','latex');
% ylabel('$rv$/rad', 'interpreter','latex');
% % yticks([0,.2, .4, .6, .8, 1.0, 1.2, 1.4]);
% set(gca,'FontSize', 36);
% 
figure;
plot(t, pos_actual, t, pos_d, '--', 'LineWidth', 2);
xlabel("$t$/s", 'interpreter','latex');
ylabel('$p_2$/m', 'interpreter','latex');
% yticks([0,.2, .4, .6, .8, 1.0, 1.2, 1.4]);
set(gca,'FontSize', 36);
% 
figure;
plot3(pos_actual(1,:),pos_actual(2,:),pos_actual(3,:), 'r-', pos_d(1,:),pos_d(2,:),pos_d(3,:), 'b--', 'LineWidth', 2);
xlabel("$x$/m", 'interpreter','latex');
ylabel('$y$/m', 'interpreter','latex');
zlabel("$z$/m", 'interpreter','latex');
% xticks([0.5, 0.54, 0.58]);
set(gca,'FontSize', 36);
axis equal;
grid on;
end

function [u, e] = speed_controller(t, y, desired_speed, Kp, Ki)
% y(1): theta, y(2) : w, y(3) sum e
n = length(y) / 4;
[wd, e2] = desired_speed(t, y);
e1 = wd - y(n + 1: 2 * n);
u = Kp .* e1 + Ki .* y(2 * n + 1 : 3 * n);
e = [e1;e2];
end

function [desired_speed, e] = task_rcm_controller(t, y, desired_rcm_pos, robot, Ts, p1, p2, kesai, Kp, Ki, Kd)
Rs = Ts(1:3,1:3);
ps = Ts(1:3,4);
Tp1 = [eye(3), -[0, 0, p1]'; 0 0 0 1];
Tp2 = [eye(3), -[0, 0, p2]'; 0 0 0 1];

n = (length(y)) / 4;
q = y(1:n);
% qd = y(n + 1 : 2*n);
lambda = y(end);
Jb = jacobian_matrix(robot, q);
T = forward_kin_general(robot, q);
q0 = inverse_kin_kuka_robot_kesai_near(robot, T, kesai, q)';
R = T(1:3,1:3);
p = T(1:3,4);
J1 = adjoint_T(Tp1) * Jb;
J2 = adjoint_T(Tp2) * Jb;
J1 = R * J1(4:6,:);
J2 = R * J2(4:6,:);
P1 = p + R * [0, 0, p1]';
P2 = p + R * [0, 0, p2]';
tipxyz = Rs' * (P2 - ps);
rcm = P1 + (P2 - P1) * lambda;
Jrcm = [J1 + lambda * (J2 - J1), P2 - P1];
Jt = Rs' * J2;
[xyz, prcm] = desired_rcm_pos(t);
e = [xyz - tipxyz;prcm - rcm];
J = [Jt, zeros(3,1); Jrcm];
td = Kp(1:3) .* e(1:3) + Ki(1:3) .* y(3*n + 1 : 3*n + 3);
rcmd = Kp(4:6) .* e(4:6) + Ki(4:6) .* y(3*n + 4 : 3*n + 6);
% 零空间控制
w = [q0 - q;0];
desired_speed = lsqminnorm(J,[[td;rcmd],J * w]);
desired_speed = desired_speed(:,1) + w - desired_speed(:,2);
e = [e;desired_speed(end)];
desired_speed = desired_speed(1:n);
end



function [xyz, prcm]= desired_rcm_pos(t, Ts, p1, p2, lambda0)
%相对于初始位置 
xyz = [0;0.02*t; p2];
% xyz = [0.02*sin(2*t);0.02*t; p2];
% xyz = [0.03*cos(t) - 0.03; 0.03 * sin(t); p2];
% xyz = [0.02*cos(t); 0.02 * sin(t); p2 + 0.005*t];
% xyzd = [-0.02*sin(t); 0.02 * cos(t); 0];
% xyzd = zeros(3,1);
Rs = Ts(1:3,1:3);
ps = Ts(1:3,4);
P1 = ps + Rs * [0, 0, p1]';
P2 = ps + Rs * [0, 0, p2]';
prcm = P1 + (P2 - P1) * lambda0;% + [0,0,0.01]'; % change rcm
end

function ret = odeplot(t, y, flag, port)
if strcmp(flag, 'init') == 1
elseif isempty(flag)
    n = size(y,1) / 4;
    % n = length(y) / 4;
    setJoints(port, y(1:n, end));
    % torque = [torque, y(2 * n+1:3*n, end)];
else
end
ret = 0;
end

function setJoints(port, jt)
cmd = sprintf('robot;%f;%f;%f;%f;%f;%f;%f;', jt(1), jt(2), jt(3), jt(4), jt(5)...
    ,jt(6), jt(7));
writeline(port,cmd,"127.0.0.1",7755);
end

function ptp(port, jts, vel)
if nargin < 3
    vel = .4;
end
start = queryJoints(port);
wayPoints = [start',jts'];
Freq = 200;
rate = rateControl(Freq);
T = max(abs(jts - start) / vel);
numSamples = round(T * Freq) + 1;
jt = trapveltraj(wayPoints,numSamples);
for i = 1 : numSamples
    setJoints(port, jt(:,i));
    waitfor(rate);
end
end

function joints = queryJoints(port)
% ";" 表示查询关节角
writeline(port,"robot;","127.0.0.1",7755);
s = readline(port);
joints = sscanf(s,'%f;%f;%f;%f;%f;%f;%f;')';
joints = mod(joints + pi, 2*pi) - pi;
end

function yd = joint_motor_dynamic(t, y, u, d, J, B, r)
% y(1): theta, y(2): w, y(3): sum e
n = length(y) / 4;
yd = zeros(4 * n, 1);
[U, e] = u(t,y);
yd(1:n) = y(n +1 : 2 * n);
yd(n +1 : 2*n) = (U - d(t, y) ./ r - B .* y(n +1: 2 *n)) ./ J;
yd(2 * n + 1 : end) = e;
% disp(u(t,y));
end

function dist = dist2line(pt, p1, p2)
    s = p2 - p1;
    s = s / norm(s);
    l1 = pt - p1;
    l2 = dot(l1, s);
    dist = real(sqrt(l1'*l1 - l2*l2));
end
