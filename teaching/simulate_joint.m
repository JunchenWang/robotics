function simulate_joint
joint = joint_model(0.1, 0.2, 200);
y0 = [0, 0]';
u = @(t, y) joint.B * 5;
d = @(t, y) 10;
dynamic = @(t, y) joint_dynamics(joint, t, y, u, d);
tspan = [0, 10];
[t, y] = ode45(dynamic, tspan, y0);
plot(t, y(:,2), 'b-','LineWidth',2);
title('关节速度开环控制 u= 5B, d = 10');
xlabel('$t$ / s','interpreter','latex');
ylabel('$\dot{q}$ / (rad/s)', 'interpreter','latex');
set(gca,'FontSize', 16);