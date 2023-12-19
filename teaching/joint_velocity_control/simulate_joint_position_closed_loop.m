function simulate_joint_position_closed_loop
joint = joint_model(0.1, 0.2, 200);
y0 = [0.5, 1, 0, 0]';
Kp = 1;
Ki = 2;
Kp_pos = 10;
Ki_pos = 0;
Kd_pos = 1;
p = @(t, y) pos_controller(t, y, @desired_pos, Kp_pos, Ki_pos, Kd_pos);
u = @(t, y) speed_controller(t,y, p, Kp, Ki);
d = @(t, y) 10;
dynamic = @(t, y) joint_dynamics(joint, t, y, u, d);
tspan = [0, 10];
[t, y] = ode45(dynamic, tspan, y0);
q_d = desired_pos(t);
plot(t, y(:,1), 'b-', t, q_d, 'r-', 'LineWidth',2);
title('关节位置控制Kp=10,Ki=0,Kd=1');
xlabel('$t$ / s','interpreter','latex');
xticks(linspace(0,10, 11));
ylabel('$q$ / (rad)', 'interpreter','latex');
set(gca,'FontSize', 16);
lg = legend('仿真位置','期望位置');
fontsize(lg,14,'points')
end
function [qd, qe] = pos_controller(t, y, q_d, Kp, Ki, Kd)
    [q_d, qd_d] = q_d(t);
    qe = q_d - y(1);
    qed = qd_d - y(2);
    sum_qe = y(4);
    qd = qd_d + Kp * qe + Ki * sum_qe + Kd * qed; 
end
function [u, qed, qe] = speed_controller(t, y, qd_d, Kp, Ki)
    qd = y(2);
    [qd_d, qe] = qd_d(t, y);
    qed = qd_d - qd;
    sum_qed = y(3);
    u = Kp * qed + Ki * sum_qed;
end
function [q_d, qd_d] = desired_pos(t)
    q_d = 2 * sin(t);
    qd_d = 2 * cos(t);
    % q_d = 2 * ones(size(t));
    % qd_d = 0;
end
function yd = joint_dynamics(joint, t, y, u, d)
% y: q, qd, qed的积分, qe的积分
yd = zeros(4,1);
[u, qed, qe] = u(t,y);
yd(1) = y(2);
yd(2) = (u - d(t, y) / joint.r - joint.B * y(2)) / joint.J;
yd(3) = qed;
yd(4) = qe;
end
