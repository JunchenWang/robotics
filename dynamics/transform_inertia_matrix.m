function I = transform_inertia_matrix(Ib, m, Tb)
% Ib defined in a center mass frame
% Tb is with respect to I's frame
R = Tb(1:3,1:3);
t = Tb(1:3,4);
I = R*Ib*R' + m * (t'*t *eye(3) - t * t');