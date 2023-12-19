function simulate_two_link_robot_joint_space
port = udpport("byte");
robot = two_link_robot;
n = robot.dof;
y0 = zeros(n * 2, 1);
y0(1:n) = y0(1:n) + rand(2,1) * 1e-3;
ptp_move(port, y0(1:n));
Kp = eye(2);
Kd = 2*eye(2);
tau = @(t, y) joint_controller(t, y, @desired_pos, robot, Kp, Kd);
tspan = [0, 10];
opts = odeset('OutputFcn',@(t, y, flag) odeplot(t, y, flag, robot, port));
dynamic = @(t, y) two_link_dynamics(robot, t, y, tau);
[t, y] = ode45(dynamic, tspan, y0, opts);
len = length(t);
q_d = zeros(len, 2);
for i = 1 : len
    q_d(i,:) = desired_pos(t);
end
plot(t, y(:,1), 'b-*', t, y(:,2), 'b-s', t, q_d, 'r-', 'LineWidth',2);
title('双连杆机器人关节空间控制');
xlabel('$t$ / s','interpreter','latex');
xticks(linspace(0,10, 11));
ylabel('$q$ / (rad)', 'interpreter','latex');
set(gca,'FontSize', 16);
lg = legend('仿真位置 q1','仿真位置 q2', '期望位置 q1', '期望位置 q2');
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
function tau = joint_controller(t, y, q_d, robot, Kp, Kd)
q = y(1:2);
qd = y(3:4);
[q_d, qd_d, qdd_d] = q_d(t);
qe = q_d - q;
qed = qd_d - qd;
M = mass_matrix(robot, q);
C = c_matrix(robot, q, qd);
G = g_matrix(robot, q);
damping = robot.b * qd;
tau = M * (qdd_d + Kp * qe + Kd * qed) + C * qd + G + damping;
end
function [q_d, qd_d, qdd_d] = desired_pos(t)
    q_d = [pi / 4, pi / 3]';
    qd_d = zeros(2,1);
    qdd_d = zeros(2,1);
end