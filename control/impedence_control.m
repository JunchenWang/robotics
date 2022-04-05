function [re, pe] = impedence_control(robot, Xd, Vd, Mp, Bp, Kp, Mr, Br, Kr, y, f, dt, re, pe)
n = robot.dof;
q = y(1:n);
qd = y(n + 1 : 2 * n);
% velocity noise
% qd = max(qd) * randn(n,1); 
% X = forward_kin_kuka(q);
% X(1:3,4) = X(1:3,4) / 1000;
[Jb, X] = jacobian_matrix(robot, q);
% disp(norm(X-XX));
R = X(1:3,1:3);
% p = X(1:3,4);
Rd = Xd(1:3, 1:3);
% pd = Xd(1:3,4);
V = Jb * qd;
% pe = R' * (pd - p);
ped = -so_w(V(1:3)) * pe + R' * Vd(4:6) - V(4:6);
% re = logR(R'*Rd)';
Ar = w_dr_A(re);
red = pinv(Ar) * (Rd' * Vd(1:3) - Rd'*R*V(1:3));
pedd = Mp \ (f(4:6) - Bp * ped - Kp * pe);
redd = Mr \ (f(1:3) - Br * red - Kr * re);

% pe = pe + ped * dt + 0.5 * pedd * dt^2;
ped = ped + dt * pedd;
pe = pe + dt * ped;

% re = re + red * dt + 0.5 * redd * dt^2;
red = red + dt * redd;
re = re + dt * red;

% R = Rd * exp_w(re)';
% p = pd - R * pe;
% 
% X = [R, p; 0, 0, 0, 1];

