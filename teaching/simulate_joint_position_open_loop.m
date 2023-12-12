function simulate_joint_position_open_loop
joint = joint_model(0.1, 0.2, 200);
y0 = [0.5, 0.1, 0]';
Kp = 1;
Ki = 2;
u = @(t, y) speed_controller(t,y, @desired_velocity, Kp, Ki);
d = @(t, y) 10;
dynamic = @(t, y) joint_dynamics(joint, t, y, u, d);
tspan = [0, 10];
[t, y] = ode45(dynamic, tspan, y0);
q_d = desired_pos(t);
plot(t, y(:,1), 'b-', t, q_d, 'r-', 'LineWidth',2);
title('关节位置控制——开环');
xlabel('$t$ / s','interpreter','latex');
xticks(linspace(0,10, 11));
ylabel('$q$ / (rad)', 'interpreter','latex');
set(gca,'FontSize', 16);
lg = legend('仿真位置','期望位置');
fontsize(lg,14,'points')
end

function [u, qed] = speed_controller(t, y, qd_d, Kp, Ki)
    qd = y(2);
    qed = qd_d(t) - qd;
    sum_qed = y(3);
    u = Kp * qed + Ki * sum_qed;
end
function [q_d, qd_d] = desired_pos(t)
    q_d = 2 * sin(t);
    qd_d = 2 * cos(t);
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
    % qd_d = 5 * ones(size(t));
    [~,qd_d] = desired_pos(t);
end