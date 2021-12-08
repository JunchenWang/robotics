function yd = motorDynamics(t, y, u, motor)
% y(1): angle, y(2): angular vel, y(3): angular acc
R = motor.R;
L = motor.L;
b = motor.b;
K = motor.K;
J = motor.J;
yd = zeros(3,1);
yd(1) = y(2);
yd(2) = y(3);
yd(3) = (u(t, y)-(R*J + L*b)/K*y(3) - (R*b/K + K)*y(2)) * K/(L *J);
% disp(u(t,y));