function simulate_joint_closed_loop
joint = joint_model(0.1, 0.2, 200);
y0 = [0, 0, 0]';
Kp = 2;
Ki = 4;
u = @(t, y) speed_controller(t,y, @desired_velocity, Kp, Ki);
d = @(t, y) 10*sin(t);
dynamic = @(t, y) joint_dynamics(joint, t, y, u, d);
tspan = [0, 10];
[t, y] = ode45(dynamic, tspan, y0);
plot(t, y(:,2), 'b-', t, desired_velocity(t), 'r--', 'LineWidth',2);
title('关节速度闭环控制');
xlabel('$t$ / s','interpreter','latex');
xticks(linspace(0,10, 11));
ylabel('$\dot{q}$ / (rad/s)', 'interpreter','latex');
% yticks(linspace(0,6, 7));
% ylim([0,5.5]);
set(gca,'FontSize', 32);
lg = legend('仿真速度','期望速度');
fontsize(lg,32,'points');
set(gcf,'Position',[100 100 1200 800]);
end

function [u, qed] = speed_controller(t, y, qd_d, Kp, Ki)
    qd = y(2);
    qed = qd_d(t) - qd;
    sum_qed = y(3);
    u = Kp * qed + Ki * sum_qed;
end

function yd = joint_dynamics(joint, t, y, u, d)
% y: q, qd, qed的积分
yd = zeros(3,1);
[u, e] = u(t,y);
yd(1) = y(2);
yd(2) = (u - d(t, y) / joint.r - joint.B * y(2)) / joint.J;
yd(3) = e;
end

function qd_d = desired_velocity(t)
    qd_d = 5 * ones(size(t));
    % qd_d = sin(t);
end