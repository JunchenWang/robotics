function simulate_two_link_robot_task_space
port = udpport("byte");
robot = two_link_robot;
n = robot.dof;
y0 = zeros(n * 2, 1);
y0(1:n) = two_link_IK(robot, [0.6, 0], 1);
ptp_move(port, y0(1:n));
tspan = [0, 10];
Kp = 10*eye(2);
Kd = 1*eye(2);
[~,~,~,~,pp] = trapveltraj([0, 1], 100, 'EndTime', tspan(2));
pp = pp{1};
desired_task_pos = @(t) desired_task_pos_circle(t, pp, fnder(pp, 1), fnder(pp, 2));
tau = @(t, y) task_controller(t, y, desired_task_pos, robot, Kp, Kd);

opts = odeset('OutputFcn',@(t, y, flag) odeplot(t, y, flag, robot, port));
dynamic = @(t, y) two_link_dynamics(robot, t, y, tau);
[t, y] = ode45(dynamic, tspan, y0, opts);

len = length(t);
x_d = zeros(len, 2);
x = zeros(len, 2);
for i = 1 : len
    x_d(i,:) = desired_task_pos(t(i));
    x(i,:) = two_link_FK(robot, y(i,1:2));
end
plot(t, x(:,1), 'b-*', t, x(:,2), 'b-s', t, x_d, 'r-', 'LineWidth',2);
title('双连杆机器人任务空间控制');
xlabel('$t$ / s','interpreter','latex');
xticks(linspace(0,10, 11));
ylabel('$x-y$ / (m)', 'interpreter','latex');
set(gca,'FontSize', 16);
lg = legend('仿真位置 x','仿真位置 y', '期望位置 x', '期望位置 y');
fontsize(lg,14,'points')
end

function yd = two_link_dynamics(robot, t, y, tau)
% q, qd 是状态
q = y(1:2);
qd = y(3:4);
yd = zeros(4,1);
yd(1:2) = qd;
M = mass_matrix(robot, q);
C = c_matrix(robot, q, qd);
G = g_matrix(robot, q);
damping = robot.b * qd;
yd(3:4) = M \ (tau(t,y) - C * qd - G - damping);
end

function ret = odeplot(t, y, flag, robot, port)
if strcmp(flag, 'init') == 1
elseif isempty(flag)
    set_joints(port, y(1:robot.dof, end));
else
end
ret = 0;
end

function tau = task_controller(t, y, x_d, robot, Kp, Kd)
q = y(1:2);
qd = y(3:4);
[x_d, xd_d, xdd_d] = x_d(t);
x = two_link_FK(robot, q);
[J, dJ] = j_matrix(robot, q, qd);
xd = J * qd;
xe = x_d - x;
xed = xd_d - xd;
xdd_c = xdd_d + Kp * xe + Kd * xed;
qdd_c = J \ (xdd_c - dJ * qd);
M = mass_matrix(robot, q);
C = c_matrix(robot, q, qd);
G = g_matrix(robot, q);
damping = robot.b * qd;
tau = M * qdd_c + C * qd + G + damping;
end

function [x_d, xd_d, xdd_d] = desired_task_pos_circle(t, pp, ppd, ppdd)
    % x_d = [0, 0.6]';
    % xd_d = zeros(2,1);
    % xdd_d = zeros(2,1);
    s = ppval(pp, t);
    sd = ppval(ppd, t);
    sdd = ppval(ppdd, t);
    w = 4*pi;
    x_d = [0.6*cos(w*s), 0.6*sin(w*s)]';
    xd_d = [-0.6*w*sin(w*s)*sd, 0.6*w*cos(w*s)*sd]';
    xdd_d = [-0.6*w^2*cos(w*s)*sd^2 - 0.6*w*sin(w*s)*sdd , -0.6*w^2*sin(w*s)*sd^2 + 0.6*w*cos(w*s)*sdd]';
end

function [x_d, xd_d, xdd_d] = desired_task_pos_line(t, pp, ppd, ppdd)
    % x_d = [0, 0.6]';
    % xd_d = zeros(2,1);
    % xdd_d = zeros(2,1);
    s = ppval(pp, t);
    sd = ppval(ppd, t);
    sdd = ppval(ppdd, t);
    x_d = [0.6, 0.6 * s]';
    xd_d = [0, 0.6*sd]';
    xdd_d = [0, 0.6*sdd]';
end