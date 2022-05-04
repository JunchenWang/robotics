function [T, re, pe, red, ped, redd, pedd] = impedance_control(robot, Tcp, Tsensor, Xd, Vd, Mp, Bp, Kp, Mr, Br, Kr, y, f, dt)
n = robot.dof;
q = y(1:n);
qd = y(n + 1 : 2 * n);
[Jb, X] = jacobian_matrix(robot, q);
Jb = adjoint_T(InvertT(Tcp)) * Jb;
X = X * Tcp;
R = X(1:3,1:3); 
p = X(1:3,4);
Rd = Xd(1:3, 1:3);
pd = Xd(1:3,4);
V = Jb * qd;
pe = R' * (pd - p);
ped = -so_w(V(1:3)) * pe + R' * Vd(4:6) - V(4:6);
[Um, ~, Vm] = svd(R'*Rd);
tem = Um * Vm';
re = logR(tem)';
% disp([norm(pe2 - pe) / norm(pe), norm(re2 - re) / norm(re)]);
Ar = w_dr_A(re);
f = adjoint_T(InvertT(Tsensor)*Tcp)'*f;
red = Ar \ (Rd' * Vd(1:3) - Rd'*R*V(1:3));

pedd = Mp \ (f(4:6) - Bp * ped - Kp * pe);
redd = Mr \ (f(1:3) - Br * red - Kr * re);
% pe = pe + ped * dt + 0.5 * pedd * dt^2;
ped = ped + dt * pedd;
pe = pe + dt * ped;
% re = re + red * dt + 0.5 * redd * dt^2;
red = red + dt * redd;
re = re + dt * red;
Tbs = eye(4);
Tbs(1:3,1:3) = Xd(1:3,1:3)';
Vd_b = adjoint_T(Tbs) * Vd * dt;
Td = Xd * exp_twist(Vd_b');
R = Td(1:3,1:3) * exp_w(-re');
p = Td(1:3,4) - R * pe;
T = [R, p; 0, 0, 0, 1];


