function tao = setpoint_controller(robot, t, y, Td, veld, accd, Kp, Kd)
% error in desired frame
n = robot.dof;
q = y(1:n);
qd = y(n + 1 : 2 * n);
[M, C, G, Jb, dJb, T] = mass_c_g_matrix(robot, q, qd);
Vb = Jb * qd;
wb = Vb(1:3);
vb = Vb(4:6);
wd = veld(1:3);
dwd = accd(1:3);
vd = veld(4:6);
dvd = accd(4:6);
R = T(1:3,1:3);
p = T(1:3, 4);
Rd = Td(1:3, 1:3);
pd = Td(1:3, 4);
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

% invJ = pinv(Jx);
% Vt = [invAr * R' * wd; Rd'*vd - so_w(pe) * Rd' * wd];

dJx = [dinvAr * Jw + invAr * dJw; dRd' * R * Jv + Rd' * dR * Jv + Rd' * R * dJv];
dVt = [dinvAr * R' * wd + invAr * dR' * wd + invAr * R' * dwd;
       dRd' * vd + Rd'*dvd - so_w(ped) * Rd' * wd - so_w(pe) * dRd' * wd - so_w(pe) * Rd' * dwd];
xe = [re;pe];
xed = [red; ped];
ax = dVt - Kp * xe - Kd * xed;
aq = Jx \ (ax - dJx * qd);
tao = G + M * aq + C * qd;
% disp(norm(dJx - dJb));
