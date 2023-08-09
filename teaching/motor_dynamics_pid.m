function yd = motor_dynamics_pid(t, y, u)
% 系统第三个状态量是累计误差，它的导数是当前误差e，由控制器返回
R = 0.6;
L = 0.12;
b = 1e-1;
K = 25.9e-3;
J = 0.1;
yd = zeros(3,1);
yd(1) = y(2);
[U, e] = u(t, y);
yd(2) = (U-(R*J + L*b)/K*y(2) - (R*b/K + K)*y(1)) * K/(L *J);
yd(3) = e;