function [u, e] = pid_controller(t, y, yd, Kp, Ki, Kd)
ydt = yd(t);
e = ydt(1) - y(1);
ed = ydt(2) - y(2);
acc_e = y(3);
u = Kp * e + Kd * ed + Ki * acc_e;

