function simulate_motor
y0 = [0, 0]';
u = @(t, y) 12;
dynamic = @(t, y) motor_dynamics(t, y, u);
tspan = [0, 10];
[t, y] = ode45(dynamic, tspan, y0);
plot(t, y');