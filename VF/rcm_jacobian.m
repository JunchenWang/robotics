function [J, dJ, x] = rcm_jacobian(robot, q, dq, p1_F, p2_F, rcm)
% x: the nearest point on p1-p2 to rcm. rcm is defined in base frame.
T1 = [eye(3), p1_F; 0 0 0 1];
T2 = [eye(3), p2_F; 0 0 0 1];
[dJb, Jb, dT, T] = derivative_jacobian_matrix(robot, q, dq);
R = T(1:3, 1:3);
dR = dT(1:3,1:3);
t = T(1:3,4);
J1_all = adjoint_T(tform_inv(T1)) * Jb;
dJ1_all =  adjoint_T(tform_inv(T1)) * dJb;
J1 = R * J1_all(4:6,:);
dJ1 = dR * J1_all(4:6,:) + R * dJ1_all(4:6,:);

J2_all = adjoint_T(tform_inv(T2)) * Jb;
dJ2_all =  adjoint_T(tform_inv(T2)) * dJb;
J2 = R * J2_all(4:6,:);
dJ2 = dR * J2_all(4:6,:) + R * dJ2_all(4:6,:);

P1 = R * p1_F + t;
P2 = R * p2_F + t;
v = P2 - P1;
u = rcm - P1;
dv = (J2 - J1)*dq;
du = -J1 * dq;
vn2 = dot(v,v);
vn4 = vn2 * vn2;
lambda = dot(u,v) / vn2;
x = P1 + lambda * (P2 - P1);
term1 = 2*dot(u,v)*v' - (u + v)' / vn2;
term2 = u' / vn2 - 2 * dot(u,v) * v';
J_lambda = term1 * J1 + term2 * J2;
dterm1 = 2 * (dot(du,v) + dot(u,dv)) * v' + 2 * dot(u,v) * dv' - ((du + dv)' * vn2 - 2 * dot(v, dv) * (u + v)') / vn4;
dterm2 = (du' * vn2 - 2 * dot(v, dv) * u') / vn4 - (2 * (dot(du,v) + dot(u,dv)) * v' + 2 * dot(u,v) * dv');
dJ_lambda = dterm1 * J1 + term1 * dJ1 + dterm2 * J2 + term2 * dJ2;

J = J1 + lambda * (J2 - J1) + v * J_lambda;
dJ = dJ1 + J_lambda * dq * (J2 - J1) + lambda * (dJ2 - dJ1) + dv * J_lambda + v * dJ_lambda;







