function [T, re, pe, red, ped] = impedance_k_b_control(robot, Tcp, Tsensor, Xd, Vd, Bp, Kp, Br, Kr, q, f, dt)
X = forward_kin_general(robot, q);
X = X * Tcp;
R = X(1:3,1:3); 
p = X(1:3,4);
Rd = Xd(1:3, 1:3);
pd = Xd(1:3,4);
pe = R' * (pd - p);
re = logR(R'*Rd)';
f = adjoint_T(InvertT(Tsensor)*Tcp)'*f;
ped = Bp \ (f(4:6) - Kp * pe);
red = Br \ (f(1:3) - Kr * re);
pe = pe + dt * ped;
re = re + dt * red;
Tbs = eye(4);
Tbs(1:3,1:3) = Xd(1:3,1:3)';
Vd_b = adjoint_T(Tbs) * Vd * dt;
Td = Xd * exp_twist(Vd_b');
R = Td(1:3,1:3) * exp_w(-re');
p = Td(1:3,4) - R * pe;
T = [R, p; 0, 0, 0, 1];


