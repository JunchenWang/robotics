function Ib = transform_com_inertia_matrix(I, m, Tb)
% Tb is a center mass frame with respect to I's frame
T = tform_inv(Tb);
R = T(1:3,1:3);
t = T(1:3,4);
q = R'*t;
I_rot = I - m * (q'*q *eye(3) - q * q');
Ib = R * I_rot * R';