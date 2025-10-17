function simulate_joint_torque_control
joint = joint_dynamic_model(0.1, 0.5, 0.1, 0.2);
y0 = [0, 0, 0]';
Kp = 20;
Ki = 0;
Kd = 10;
p = @(t, y) controller(joint, t, y, @desired_pos, Kp, Ki, Kd);
dynamic = @(t, y) joint_dynamics(joint, t, y, p);
tspan = [0, 10];
[t, y] = ode45(dynamic, tspan, y0);
q_d = desired_pos(t);
plot(t, y(:,1), 'b-', t, q_d, 'r-', 'LineWidth',2);
title('关节位置控制:力矩模式');
xlabel('$t$ / s','interpreter','latex');
xticks(linspace(0,10, 11));
ylabel('$q$ / (rad)', 'interpreter','latex');
set(gca,'FontSize', 16);
lg = legend('仿真位置','期望位置');
fontsize(lg,14,'points')
end
function [tau, qe] = controller(joint, t, y, q_d, Kp, Ki, Kd)
    [q_d, qd_d, qdd_d] = q_d(t);
    qe = q_d - y(1);
    qed = qd_d - y(2);
    sum_qe = y(3);
    % tau = Kp * qe + Ki * sum_qe + Kd * qed;
    tau = joint.M * (qdd_d + Kp * qe + Ki * sum_qe + Kd * qed)...
           + joint.m * 9.8 * joint.r * cos(y(1)) + joint.b * y(2);
end

function [q_d, qd_d, qdd_d] = desired_pos(t)
    % q_d = 2 * sin(t);
    % qd_d = 2 * cos(t);
    q_d = 2 * ones(size(t));
    qd_d = 0;
    qdd_d = 0;
end
function yd = joint_dynamics(joint, t, y, tau)
% y: q, qd, qe的积分
yd = zeros(3,1);
g = 9.8;
[tau, qe] = tau(t,y);
yd(1) = y(2);
yd(2) = (tau - joint.m*g*joint.r*cos(y(1)) - joint.b * y(2)) / joint.M;
yd(3) = qe;
end
