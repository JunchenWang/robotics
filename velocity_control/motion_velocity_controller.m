function qv = motion_velocity_controller(robot, t, y, Td, veld, accd, Kp, Ki, Kd, dt)
% error in desired frame
persistent sum;
persistent pre_t;
persistent pre_xe;
if isempty(sum)
    sum = zeros(6,1);
    pre_t = t;
    pre_xe = zeros(6,1);
end
n = robot.dof;
q = y(1:n);
qd = y(n + 1 : 2 * n);
[Jb, T] = jacobian_matrix(robot, q);
Vb = Jb * qd;
wb = Vb(1:3);
vb = Vb(4:6);
wd = veld(1:3);
vd = veld(4:6);
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
Jw = Jb(1:3,:);
Jv = Jb(4:6,:);
Jx = [Ar \ Jw; Rd'*R*Jv];
Vt = [Ar \  R' * wd; Rd'*vd - so_w(pe) * Rd' * wd];

xe = [re;pe];
xed = [red; ped];
sum = sum + (t - pre_t) * (xe + pre_xe) / 2;
pre_t = t;
pre_xe = xe;
qv = pinv(Jx) * (Vt - Kp * xe - Ki * sum - Kd * xed);

