function control_motor_velocity
motor.R = 0.6;
motor.L = 0.12;
motor.b = 1e-1;
motor.K = 25.9e-3;
motor.J = 0.1;
Kp = 6;
Ki = 4;
Kd = 1;
yd = @(t) [5; 0]; %速度控制，不用指定位置
clear pidVelocityController;
u = @(t, y) pidVelocityController(t, y, yd, Kp, Ki, Kd);
y0 = [0, 0, 0]';
dynamic = @(t, y) motorDynamics(t, y, u, motor);
tspan = [0, 10];
[t, y] = ode45(dynamic, tspan, y0);
for i = 1 : length(t)
    ydt = yd(t(i));
    y(i,3) = ydt(1) - y(i,2);
end
plot(t, y(:,2:3)');