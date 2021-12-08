function simulate_motor
motor.R = 0.6;
motor.L = 0.12;
motor.b = 1e-1;
motor.K = 25.9e-3;
motor.J = 0.1;
y0 = [0, 0, 0]';
u = @(t, y) 12;
dynamic = @(t, y) motorDynamics(t, y, u, motor);
tspan = [0, 10];
[t, y] = ode45(dynamic, tspan, y0);
plot(t, y(:,2:3)');