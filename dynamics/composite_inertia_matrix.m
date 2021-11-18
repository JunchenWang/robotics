function [I, T] = composite_inertia_matrix(I1, m1, I2, m2, T12)
% combine I1 I2 to T, T is also relative to T1
t = m2 / (m1 + m2) * T12(1:3,4);
T = [eye(3), t; 0 0 0 1];
I = transform_inertia_matrix(I1, m1, T) + transform_inertia_matrix(I2, m2, tform_inv(T12)*T);
