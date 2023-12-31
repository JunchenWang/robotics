function simulate_joint_open_loop
joint = joint_model(0.1, 0.2, 200);
y0 = [0, 0]';
u = @(t) open_loop_controller(t, @desired_velocity, joint);
d = @(t) 10;
dynamic = @(t, y) joint_dynamics(joint, t, y, u, d);
tspan = [0, 10];
[t, y] = ode45(dynamic, tspan, y0);
qd_d = desired_velocity(t);
plot(t, y(:,2), 'b', t, qd_d, 'r--', 'LineWidth',2);
title('关节速度开环控制');
xlabel('$t$ / s','interpreter','latex');
xticks(linspace(0,10, 11));
ylabel('$\dot{q}$ / (rad/s)', 'interpreter','latex');    
yticks(linspace(0,6, 7));
ylim([0,5.5]);
set(gca,'FontSize', 32);
lg = legend('仿真速度','期望速度');
fontsize(lg,32,'points');
set(gcf,'Position',[100 100 1200 800]);
end

function u = open_loop_controller(t, qd_d, joint)
    [qd_d, qdd_d] = qd_d(t);
    u = joint.J * qdd_d + joint.B * qd_d;
end

function yd = joint_dynamics(joint, t, y, u, d)
% y: q, qd
yd = zeros(2,1);
yd(1) = y(2);
yd(2) = (u(t) - d(t) / joint.r - joint.B * y(2)) / joint.J;
end

function [qd_d, qdd_d] = desired_velocity(t)
    qd_d = 5 * ones(size(t));
    % qd_d = sin(t);
    qdd_d = 0;
end
