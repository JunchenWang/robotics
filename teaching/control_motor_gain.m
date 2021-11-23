function control_motor_gain
Kp = 20;
wd = 5;
u = @(t, y) Kp * (wd - y(1));
y0 = [0, 0]';
dynamic = @(t, y) motor_dynamics(t, y, u);
tspan = [0, 10];
[t, y] = ode45(dynamic, tspan, y0);
y(:,2) = wd - y(:,1);
plot(t, y');