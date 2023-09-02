function [I, T] = composite_inertias(I1, m1, I2, m2, T1, T2)
% combine I1 I2 to T, T is also relative to T1
ro = m2 / (m1 + m2);
t = ro * T2(1:3,4) + (1 - ro) * T1(1:3,4);
T = [eye(3), t; 0 0 0 1];
I = transform_inertia_matrix(I1, m1, tform_inv(T1)) + transform_inertia_matrix(I2, m2, tform_inv(T2));
