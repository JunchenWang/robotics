function Ib = transform_com_inertia_matrix(I, m, Tb)
% Tb is a center mass frame with respect to I's frame
R = Tb(1:3,1:3);
t = Tb(1:3,4);
I_rot = I - m * (t'*t *eye(3) - t * t');
Ib = R' * I_rot * R;