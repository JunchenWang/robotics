function control_motor_pid
% 如果要积分，则要创建一个新的状态量，积累误差
Kp = 6;
Ki = 5;
Kd = 0.5;
yd = @(t) [5; 0];
u = @(t, y) pid_controller(t, y, yd, Kp, Ki, Kd);
y0 = [0, 0, 0]';
dynamic = @(t, y) motor_dynamics_pid(t, y, u);
tspan = [0, 10];
[t, y] = ode45(dynamic, tspan, y0);
for i = 1 : length(t)
    ydt = yd(t(i));
    y(i,2) = ydt(1) - y(i,1);
end
plot(t, y);