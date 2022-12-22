function tao = open_controller(robot, t, y, Td, veld, accd, Dd, Kd)
n = robot.dof;
q = y(1:n);
qd = y(n + 1 : 2 * n);
[dJb, Jb, dT, T] = derivative_jacobian_matrix(robot, q, qd);
invJb = pinv(Jb);
[M, C, G] = mass_c_g_matrix(robot, q, qd);
[dInvT, invT] = derivative_tform_inv(T, dT);
R = T(1:3,1:3);
p = T(1:3, 4);
Rd = Td(1:3, 1:3);
pd = Td(1:3, 4);
re = logR(R'*Rd)';
pe = R' * (p - pd);
xe = [re; pe];
X = invT * Td;

wd = Rd' * veld(1:3);
vd = Rd' * veld(4:6);
Vd = [wd;vd];


dTd = Td * se_twist(Vd);
dX = dInvT * Td + invT * dTd;
[dAdX, AdX] = derivative_adjoint_T(X, dX);
V = Jb * qd;

alphad = Rd' * accd(1:3);
ad = Rd' * accd(4:6);
dVd = [alphad;ad - cross(wd, vd)];


dxd = AdX * Vd;
xed = dxd - V;

ddxd = dAdX * Vd + AdX * dVd;


tao = G + M * invJb * dVd + (C * invJb - M * invJb * dJb * invJb) * Vd;
