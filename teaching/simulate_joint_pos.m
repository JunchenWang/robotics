function simulate_joint_pos
joint = joint_model(0.1, 0.2, 200);
y0 = [0, 0, 0]';
Kp = 1;
Ki = 10;
u = @(t, y) controller(t,y, @(t) 5, Kp, Ki);
d = @(t, y) 10;
dynamic = @(t, y) joint_dynamics(joint, t, y, u, d);
tspan = [0, 10];
[t, y] = ode45(dynamic, tspan, y0);
plot(t, y(:,2), 'b-','LineWidth',2);
title('关节速度控制 d = 10, Kp = 1, Ki = 10');
xlabel('$t$ / s','interpreter','latex');
ylabel('$\dot{q}$ / (rad/s)', 'interpreter','latex');
set(gca,'FontSize', 16);
end

function [u, qed] = controller(t, y, qd_d, Kp, Ki)
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