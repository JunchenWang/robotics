function [J, dJ, x, dx] = orientation_jacobian(Jb, dJb, T, dT)
% r: rotation vector of rotation
R = T(1:3,1:3);
dR = dT(1:3,1:3);
wb = skew_mat2vec(R'*dR)';
v = dT(1:3,4);
Jv = R * Jb(4:6,:);
dJv = dR * Jb(4:6,:) + R * dJb(4:6,:);

x = logR(T(1:3,1:3))';
Jwb = Jb(1:3,:);
dJwb = dJb(1:3,:);
Ar = w_dr_A(x);
J = Ar \ Jwb;
dx = Ar \ wb;
dAr = derivative_Ar(x, dx);
dJ = -(Ar \ dAr) * (Ar \ Jwb) + Ar \ dJwb;

J = [J;Jv(1:2,:)];
dJ = [dJ; dJv(1:2,:)];
x = [x;T(1:2,4)];
dx = [dx; v(1:2)];

end








