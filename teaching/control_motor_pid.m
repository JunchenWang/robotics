function control_motor_pid
Kp = 5;
Ki = 5;
Kd = 0;
yd = [5; 0];
clear pid_controller;
u = @(t, y) pid_controller(t, y, yd, Kp, Ki, Kd);
y0 = [0, 0]';
dynamic = @(t, y) motor_dynamics(t, y, u);
tspan = [0, 10];
[t, y] = ode45(dynamic, tspan, y0);
y(:,2) = yd(1) - y(:,1);
plot(t, y');