function control_motor_position
motor.R = 0.6;
motor.L = 0.12;
motor.b = 1e-1;
motor.K = 25.9e-3;
motor.J = 0.1;
Kp = 6;
Ki = 4;
Kd = 0.5;
K = 0.6;
yd = @(t) [5;0]; %位置控制
clear pidVelocityController;
vController = @(t, y, yd) pidVelocityController(t, y, yd, Kp, Ki, Kd);
pController = @(t, y) positionController(t, y, yd, K, vController);
y0 = [0, 0, 0]';
dynamic = @(t, y) motorDynamics(t, y, pController, motor);
tspan = [0, 10];
[t, y] = ode45(dynamic, tspan, y0);
plot(t, y(:,1)');