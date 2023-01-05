function tao = impedance_controller1(robot, t, y, Td, veld, accd, Dd, Kd, dt)
% error in actual frame
n = robot.dof;
cnt = round(t / dt) + 1;
q = y(1:n);
qd = y(n + 1 : 2 * n);
% disp(qd);
% [dJb, Jb, ~, T] = derivative_jacobian_matrix(robot, q, qd);
[M, C, G, Jb, dJb, T] = mass_c_g_matrix(robot, q, qd);
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
% [Um, ~, Vm] = svd(R'*Rd);
% tem = Um * Vm';
% re = logR(tem)';
re = logR(R'*Rd)';
pe = R' * (pd - p);
ped = R' * (vd - v) - cross(wb, pe);
Ar = w_dr_A(re);
red = Ar \ (Rd' * (wd - w));
invAr = inv(Ar);
dAr = derivative_Ar(re, red);
dinvAr = -invAr * dAr * invAr;
dRd = so_w(wd) * Rd;
dR = so_w(w) * R;
Jw = Jb(1:3,:);
Jv = Jb(4:6,:);
dJw = dJb(1:3,:);
dJv = dJb(4:6,:);
Jx = [-invAr * Rd' * R * Jw; so_w(pe) * Jw - Jv];
invJ = pinv(Jx);
Vt = [-invAr * Rd' * wd; -R'*vd];

dJx = [-dinvAr * Rd' * R * Jw - invAr * dRd' * R * Jw - invAr * Rd' * dR * Jw - invAr * Rd' * R * dJw;
      so_w(ped) * Jw + so_w(pe) * dJw - dJv];
dVt = [-dinvAr * Rd' * wd - invAr * dRd' * wd - invAr * Rd' * dwd;
        -dR' * vd - R'* dvd];


tao = G + M * invJ * (dVt) + (C * invJ - M * invJ * dJx * invJ) * (Vt)... 
       - Jx' * (Kd * [re;pe] + Dd * [red; ped]);

