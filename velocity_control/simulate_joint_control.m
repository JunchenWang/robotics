function simulate_joint_control
% 双环控制7轴机器人轨迹，状态空间4个变量：位置，速度，速度误差积分，位置误差积分
dof = 7;
y0 = zeros(4 * dof, 1);

Kp_s = [100,100,100,100,100,100,100]';
Ki_s = [50,50,50,50,100,100,100]';

Kp_p = 5 * [1, 1, 1, 1, 1, 1, 1]';
Ki_p = 0 * [1,1,1,1,1,1,1]';
Kd_p = 0;

d = @(t, y) [7;6;5;4;3;2;1]*10;
J = [8;5;4;3;2;2;2];
B = 2;
r = 200;
p = @(t, y) pos_controller(t, y, @desired_pos, Kp_p, Ki_p, Kd_p);
u = @(t, y) speed_controller(t, y, p, Kp_s, Ki_s);
dynamic = @(t, y) joint_motor_dynamic(t, y, u, d, J, B, r);
tspan = [0, 10];
[t, y] = ode45(dynamic, tspan, y0);
pos_d = zeros(dof, length(t));
vel_d = zeros(dof, length(t));
for i = 1 : length(t)
    pos_d(:, i) = desired_pos(t(i));
    vel_d(:, i) = [1;2;3;4;5;6;7];
end
figure;
plot(t, y(:,1 : dof), t, pos_d, '--', 'LineWidth', 2);
xlabel("$t$/s", 'interpreter','latex');
ylabel('$q$/rad', 'interpreter','latex');
yticks([0,.2, .4, .6, .8, 1.0, 1.2, 1.4]);
set(gca,'FontSize', 36);
% figure;
% plot(t, y(:, dof +1 : 2 * dof), t, vel_d, '--', 'LineWidth', 2);
% xlabel("t/$s$");
% ylabel('$\dot{q}$/(rad s$^{-1}$)', 'interpreter','latex');
% yticks([0,1,2,3,4,5,6,7,8]);
% set(gca,'FontSize', 36);
end

function [u, e] = speed_controller(t, y, desired_speed, Kp, Ki)
% y(1): theta, y(2) : w, y(3) sum e
n = length(y) / 4;
[wd, e2] = desired_speed(t, y);
% wd = [1,2,3,4,5,6,7]';
% e2 = zeros(7,1);
e1 = wd - y(n + 1: 2 *n);
u = Kp .* e1 + Ki .* y(2 * n + 1 : 3 * n);
e = [e1;e2];
end

function [desired_speed, e] = pos_controller(t, y, desired_pos, Kp, Ki, Kd)
n = length(y) / 4;
[posd, veld] = desired_pos(t);
e = posd - y(1:n);
ed = veld - y(n +1 : 2 * n);
desired_speed = veld + Kp .* e + Ki .* y(3 * n +1 : 4 * n) + Kd .* ed;
end

function [pos, vel] = desired_pos(t)
pos = [10 20 30 40 50 60 70]' / 180 * pi + pi / 18 * sin(t);
vel = ones(7,1) * cos(t) * pi / 18;
end

function yd = joint_motor_dynamic(t, y, u, d, J, B, r)
% y(1): theta, y(2): w, y(3): sum e
n = length(y) / 4;
yd = zeros(4 * n, 1);
[U, e] = u(t,y);
yd(1:n) = y(n +1 : 2 * n);
yd(n +1 : 2*n) = (U - d(t, y) ./ r - B .* y(n +1: 2 *n)) ./ J;
yd(2 * n + 1 : end) = e;
end