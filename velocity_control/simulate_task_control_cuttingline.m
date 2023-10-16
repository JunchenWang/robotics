function simulate_task_control
% 双环控制7轴机器人轨迹，状态空间4个变量：位置，速度，速度误差积分，位置误差积分
port = udpport("byte");

robot = convert_robot_tree2(importrobot('urdf\iiwa7\iiwa7.urdf'));
n = robot.dof;
opts = odeset('OutputFcn', @(t, y, flag) odeplot(t, y, flag, port, robot));
y0 = zeros(3 * n + 6, 1);

Kp_s = [100,100,100,100,100,100,100]';
Ki_s = [50,50,50,50,100,100,100]';

Kp_p = 5 * [1, 1, 1, 1, 1, 1]';
Ki_p = 0 * [1,1,1,1,1,1]';
Kd_p = 0 * [1,1,1,1,1,1]';

d = @(t, y) [7;6;5;4;3;2;1]*10;
J = [8;5;4;3;2;2;2];
B = 2;
r = 200;
tspan = [0, 10];

y0(1:n) = [0 75 0 -94 0 -81 0] / 180 * pi;
ptp(port, y0(1:n)');
% kesai = cal_kuka_kesai(y0);
Ts = forward_kin_general(robot, y0);
%% load path
p = @(t) desired_task_pos_cuttingline(t, path);
v = @(t, y) task_pos_controller(t, y, p, robot, Kp_p, Ki_p, Kd_p);
u = @(t, y) speed_controller(t, y, v, Kp_s, Ki_s);
dynamic = @(t, y) joint_motor_dynamic(t, y, u, d, J, B, r);

[t, y] = ode45(dynamic, tspan, y0, opts);

pos_d = zeros(6, length(t));
pos_actual = zeros(6, length(t));
rcme = zeros(1, length(t));
for i = 1 : length(t)
    T = forward_kin_general(robot, y(i,1:7));
    P1 = T * [0;0;0;1];
    P2 = T * [0;0;0.3;1];
    rcme(i) = dist2line(prcm, P1, P2);
    Td = p(t(i));
    pos_d(:, i) = Frame_T_Converter(Td);
    pos_actual(:,i) = Frame_T_Converter(T);
end

figure;
plot(t, rcme, 'LineWidth', 2);
xlabel("$t$/s", 'interpreter','latex');
ylabel('$RMS error$/m', 'interpreter','latex');
% yticks([0,.2, .4, .6, .8, 1.0, 1.2, 1.4]);
set(gca,'FontSize', 36);

figure;
plot(t, pos_actual(1:3,:), t, pos_d(1:3,:), '--', 'LineWidth', 2);
xlabel("$t$/s", 'interpreter','latex');
ylabel('$rv$/rad', 'interpreter','latex');
% yticks([0,.2, .4, .6, .8, 1.0, 1.2, 1.4]);
set(gca,'FontSize', 36);

figure;
plot(t, pos_actual(4:6,:), t, pos_d(4:6,:), '--', 'LineWidth', 2);
xlabel("$t$/s", 'interpreter','latex');
ylabel('$p$/m', 'interpreter','latex');
% yticks([0,.2, .4, .6, .8, 1.0, 1.2, 1.4]);
set(gca,'FontSize', 36);

figure;
plot3(pos_actual(4,:),pos_actual(5,:),pos_actual(6,:), 'r-', pos_d(4,:),pos_d(5,:),pos_d(6,:), 'b--', 'LineWidth', 2);
xlabel("$x$/m", 'interpreter','latex');
ylabel('$y$/m', 'interpreter','latex');
zlabel("$z$/m", 'interpreter','latex');
xticks([0.5, 0.54, 0.58]);
set(gca,'FontSize', 36);
axis equal;
grid on;
end

function [u, e] = speed_controller(t, y, desired_speed, Kp, Ki)
% y(1): theta, y(2) : w, y(3) sum e
n = (length(y) + 1) / 4;
[wd, e2] = desired_speed(t, y);
e1 = wd - y(n + 1: 2 * n);
u = Kp .* e1 + Ki .* y(2 * n + 1 : 3 * n);
e = [e1;e2];
end

function [desired_speed, e] = task_pos_controller(t, y, desired_task_pos, robot, Kp, Ki, Kd)
n = (length(y) + 1) / 4;
q = y(1:n);
qd = y(n + 1 : 2*n);
Jb = jacobian_matrix(robot, q);
V = Jb * qd;
T = forward_kin_general(robot, q);
R = T(1:3,1:3);
p = T(1:3,4);
[Td, Vd] = desired_task_pos(t);
wd = Vd(1:3);
vd = Vd(4:6);
Rd = Td(1:3,1:3);
pd = Td(1:3,4);
xe = logR(R'*Rd)';
pe = pd - p;
wb = R' * wd + Kp(1:3) .* xe + Ki(1:3) .* y(3*n + 1 : 3*n + 3) + Kd(1:3) .* (R' * wd - V(1:3));
vb = R'*(vd + Kp(4:6) .* pe + Ki(4:6) .* y(3*n + 4 : 3*n + 6) + Kd(4:6) .* (vd - R * V(4:6)));
e = [xe; pe];
desired_speed = lsqminnorm(Jb,[wb;vb]);
% disp(Jb * pinv(Jb));
end


function [Td, Vd] = desired_task_pos(t, Ts, Te, pp, ppd, ppdd)
s = ppval(pp, t);
sd = ppval(ppd, t);
sdd = ppval(ppdd, t);
ps = Ts(1:3,4);
pe = Te(1:3,4);
Rs = Ts(1:3,1:3);
Re = Te(1:3,1:3);
pd = ps + (pe - ps)*s;
xe = logR(Rs'*Re)';
Rd = Rs * exp_w(xe*s);
Td = [Rd, pd; 0 0 0 1];
Vd = [Rs * xe; (pe - ps)] * sd;
end


function [Td, Vd] = desired_task_pos_cuttingline(t, path)
% dt = 0.002;
% cnt = round(t / dt) + 1;
% Td = path()
Vd = zeros(6,1);
end




function ret = odeplot(t, y, flag, port, robot)
if strcmp(flag, 'init') == 1
elseif isempty(flag)
    setJoints(port, y(1:robot.dof, end));
else
end
ret = 0;
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

function setJoints(port, jt)
cmd = 'robot;' + join(string(jt),';') + ';';
writeline(port,cmd,"127.0.0.1",7755);
end

function joints = queryJoints(port)
writeline(port,"robot;","127.0.0.1",7755);
s = split(readline(port), ';');
joints = double(s(1:end-1))';
end

function yd = joint_motor_dynamic(t, y, u, d, J, B, r)
% y(1): theta, y(2): w, y(3): sum e
n = (length(y) + 1) / 4;
yd = zeros(3 * n + 6, 1);
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