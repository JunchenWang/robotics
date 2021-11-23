function control_motor_open
R = 0.6;
L = 0.12e-3;
b = 1e-1;
K = 25.9e-3;
J = 0.1;
wd = 5;
u = @(t, y) (R*b/K + K) * wd;
y0 = [0, 0]';
dynamic = @(t, y) motor_dynamics(t, y, u);
tspan = [0, 10];
[t, y] = ode45(dynamic, tspan, y0);
plot(t, y');