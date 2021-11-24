function yd = motor_dynamics(t, y, u)
R = 0.6;
L = 0.12;
b = 1e-1;
K = 25.9e-3;
J = 0.1;
yd = zeros(2,1);
yd(1) = y(2);
yd(2) = (u(t, y)-(R*J + L*b)/K*y(2) - (R*b/K + K)*y(1)) * K/(L *J);

