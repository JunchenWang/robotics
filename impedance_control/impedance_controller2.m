function tao = impedance_controller2(robot, t, y, Td, veld, accd, Dd, Kd, dt)
% error in desired frame
n = robot.dof;
cnt = round(t / dt) + 1;
q = y(1:n);
qd = y(n + 1 : 2 * n);
[dJb, Jb, ~, T] = derivative_jacobian_matrix(robot, q, qd);
[M, C, G] = mass_c_g_matrix(robot, q, qd);
Vb = Jb * qd;
wb = Vb(1:3);
vb = Vb(4:6);
wd = veld(1:3, cnt);
dwd = accd(1:3, cnt);
vd = veld(4:6, cnt);
dvd = accd(4:6, cnt);
R = T(1:3,1:3);
p = T(1:3, 4);
Rd = Td(1:3, 1:3, cnt);
pd = Td(1:3, 4, cnt);
w = R * wb;
v = R * vb;
re = logR(Rd'*R)';
pe = Rd' * (p - pd);
ped = Rd' * (v - vd) - cross(Rd'*wd, pe);
Ar = w_dr_A(re);
red = Ar \ (R' * (w - wd));
invAr = inv(Ar);
dAr = derivative_Ar(re, red);
dinvAr = -invAr * dAr * invAr;
dRd = so_w(wd) * Rd;
dR = so_w(w) * R;
Jw = Jb(1:3,:);
Jv = Jb(4:6,:);
dJw = dJb(1:3,:);
dJv = dJb(4:6,:);
Jx = [invAr * Jw; Rd'*R*Jv];
invJ = pinv(Jx);
Vt = [invAr * R' * wd; Rd'*vd - so_w(pe) * Rd' * wd];

dJx = [dinvAr * Jw + invAr * dJw; dRd' * R * Jv + Rd' * dR * Jv + Rd' * R * dJv];
dVt = [dinvAr * R' * wd + invAr * dR' * wd + invAr * R' * dwd;
       dRd' * vd + Rd'*dvd - so_w(ped) * Rd' * wd - so_w(pe) * dRd' * wd - so_w(pe) * Rd' * dwd];

% A = pinv(Jx') * M * invJ;
% Q = sqrtm(A);
% B = inv(Q) * Kd * inv(Q);
% Dd = 2 * Q * diag(eig(B)/10) * Q;
tao = G + M * invJ * (dVt) + (C * invJ - M * invJ * dJx * invJ) * (Vt)... 
       - Jx' * (Kd * [re;pe] + Dd * [red; ped]);

