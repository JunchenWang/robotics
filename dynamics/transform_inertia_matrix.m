function I = transform_inertia_matrix(Ib, m, T)
% Ib defined in a center mass frame
% T is with respect to com frame
R = T(1:3,1:3);
t = T(1:3,4);
I_rot = R'*Ib*R;
q = R'*t;
I = I_rot + m * (q'*q *eye(3) - q * q');