function simulate_joint_open_loop
joint = joint_model(0.1, 0.2, 200);
y0 = [0, 0]';
u = @(t) 1;
d = @(t) 10;
dynamic = @(t, y) joint_dynamics(joint, t, y, u, d);
tspan = [0, 10];
[t, y] = ode45(dynamic, tspan, y0);
plot(t, y(:,2), 'b-', 'LineWidth',2);
title('关节速度控制——开环');
xlabel('$t$ / s','interpreter','latex');
xticks(linspace(0,10, 11));
ylabel('$\dot{q}$ / (rad/s)', 'interpreter','latex');
set(gca,'FontSize', 16);
end

function yd = joint_dynamics(joint, t, y, u, d)
% y: q, qd
yd = zeros(2,1);
yd(1) = y(2);
yd(2) = (u(t) - d(t) / joint.r - joint.B * y(2)) / joint.J;
end
