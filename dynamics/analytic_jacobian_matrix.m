function Ja = analytic_jacobian_matrix(Jb, T)
R = T(1:3,1:3);
axang = rotm2axang(R);
r = axang(1:3)*axang(4);
norm_r = norm(r);
sr = vec2skew_mat(r);
Ar = eye(3) - (1 - cos(norm_r)) / norm_r^2 * sr + (norm_r - sin(norm_r))/norm_r^3 * sr*sr;
Ja = [inv(Ar), zeros(3); zeros(3), R] * Jb;